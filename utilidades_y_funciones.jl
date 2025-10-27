module Utilidades


using ForwardDiff, QuadGK
using LinearAlgebra
using SpecialPolynomials
using Polynomials
using Plots, LaTeXStrings




export chebyshevT_basis, chebyshevU_basis, matriz_G_Chebyshev,productoInt_Chebyshev, d1, d2
export Y_k_aproximado,  operador_angular, matriz_AGalerkin, sol_espectral, PL_m1, operador_radial
export  grafica5_modos_comparativa, despejar_frecuencia, R_ansatz, iteracion_newton, calculo_error_abs
export  actualizar_convergencia, grafica_convergencia

#=====================================================#
# ------------ Base  -------------

function chebyshevT_basis(N::Int)
    """funciones [T₀(x), T₁(x), ..., T_N(x)]"""
    T = Vector{Function}(undef, N+1)
    T[1] = x -> 1.0
    if N >= 1
        T[2] = x -> x
        for n in 2:N
            T[n+1] = x -> 2x*T[n](x) - T[n-1](x)
        end
    end
    return T
end

function chebyshevU_basis(N::Int)
    """funciones [U₀(x), U₁(x), ..., U_N(x)]"""
    U = Vector{Function}(undef, N+1)
    U[1] = x -> 1.0
    if N >= 1
        U[2] = x -> 2x
        for n in 2:N
            U[n+1] = x -> 2x*U[n](x) - U[n-1](x)
        end
    end
    return U
end

function matriz_G_Chebyshev(N::Int, tipo::Symbol )
    """Matriz G de  Ac=lamda Gc """
    G = (pi/2) * I(N+1)   #chebyshev U_n
    if tipo == :i   #chebyshev T_n
         G[1,1] = pi    
    end
    return G
end

function productoInt_Chebyshev(f::Function, g::Function, tipo::Symbol)
    """Matriz G de  Ac=lamda Gc """
    peso(x) = 1.0
    if tipo == :i   #chebyshev T_n
          peso = x -> 1/sqrt(1 - x^2)
          lim0=0;lim1=1
    elseif tipo==:ii  #chebyshev U_n
          peso = x -> sqrt(1 - x^2)
          lim0=-1;lim1=1
    end
   
    h(x) = g(x) * f(x) * peso(x)
    Int, _ = quadgk(h, lim0, lim1) 
    
    return Int
  
end

#=====================================================#
# ------------ Operadores -------------

function d1(y, x, h=1e-6)
    """Primera derivada dif  finitas centradas"""
    (y(x+h) - y(x-h)) / (2h)
end
function d2(y, x, h=1e-6)
    """Segunda Derivada  dif  finitas  centradas"""
    (y(x+h) - 2*y(x) + y(x-h)) / (h^2)
end


function operador_angular(w::Number, y::Function,  a::Number, s::Number, m::Number , M::Number)
    """Operador angular de la ecuación de Teukolsky"""
    c=a*w
    dy = x -> d1(y,x)
    d2y = x -> d2(y,x)
    return x -> -(1-x^2)*d2y(x)+2x*dy(x)-((c*x^2)-2*c*s+s-(m+s*x)^2/(1-x^2))*y(x)
end

function operador_radial( w::Number, H::Function,  a::Number, s::Number, m::Number , M::Number ) 
    """Operador angular  con el cambio de variable a z en [-1,1]"""
   
    d1H = z -> d1(H,z)
    d2H = z -> d2(H,z)
    i = 0 + 1im       # definimos  el numero imaginario 

    # los horizontes
    rmas= M + sqrt(M^2-a^2)
    rmenos= M - sqrt(M^2-a^2)
    R_pm= rmas-rmenos 
    
    #elemntos del cambio de  variable z=(r-r+)/(r-r-)
    r= z-> (z*rmenos-rmas)/(z-1)
    d1r= z->d1H(z)*(z-1)^2/(R_pm) 
    d2r= z->d2H(z)*(z-1)^2/(R_pm) + d1H(z)*2*(z-1)/(R_pm)
    
    #elementos del ansatz ecu.5.20
    sigmaMas=(2*w*M*rmas-m*a)/(R_pm)
    sigmaMenos=(2*w*M*rmenos-m*a)/(R_pm)
    
    # se toma el  signo positivo en el pm
    zeta=i*w
    xi=i*sigmaMas  
    eta=-i*sigmaMenos
    # se toma el signo negativo en el pm 
    #zeta=-i*w
    #xi=-s-i*sigmaMas  
    #eta=-s+i*sigmaMenos

    # elementos del operador  ecu.5.21
    Q=z-> eta/(r(z)-rmenos)+ zeta + xi/(r(z)-rmas)
    B=z-> (r(z))^2-2*M*r(z)+a^2
    K=z-> ((r(z))^2+a^2)*w-a*m
    rM=z-> r(z)-M
    S1=s+1

    termino1d= z-> B(z)*Q(z) + S1*rM(z)
    termino0d= z->2*S1*rM(z)*Q(z) + B(z)*(Q(z))^2 + ((K(z))^2 - 2*i*s*rM(z)*K(z))/(B(z)) + 4*i*s*w*r(z)

    return z-> B(z)*d2r(z) + termino1d(z)*2*d1r(z) + termino0d(z)*H(z)

end



#=====================================================#
#----------------  Metodo  Espectral -----------------

function matriz_AGalerkin(operador::Function, w::Number, tipo::Symbol, a::Number,  s::Number, m::Number , M::Number, n::Number)
    """Matriz de productos internos ⟨yi, L[yj]⟩"""
    
    if  tipo ==:i
        Basefunc = chebyshevT_basis(n-1)
    elseif tipo ==:ii
        Basefunc = chebyshevU_basis(n-1)
    end

    A = zeros(ComplexF64, n, n) 
    for (i, yi) in enumerate(Basefunc)
        for (j, yj) in enumerate(Basefunc)
            L_yj = operador(w, yj, a, s, m , M)
            Int = productoInt_Chebyshev(L_yj, yi, tipo)
            A[i, j] = Int
        end
    end
    return A
end
    

function sol_espectral(operador::Function, w::Number, tipo::Symbol, a::Number, s::Number, m::Number , M::Number, n::Int )
    """Usa el metodo  Garlekin  para  sol el problema de autofunciones"""
    
    G=matriz_G_Chebyshev(n-1, tipo)
    A=matriz_AGalerkin(operador, w, tipo, a,  s, m , M, n)
    autoval, autovec =eigen(A,G)
    return  autoval, autovec
    
end

function Y_k_aproximado(V::Matrix, tipo::Symbol, k::Int)
    """Construye la  función  resultado con los coeficientes (auto vectores) sobre la base"""

    n = size(V, 1)
    if  tipo ==:i
        Basefunc = chebyshevT_basis(n-1)
    elseif tipo ==:ii
        Basefunc = chebyshevU_basis(n-1)
    end
    return x -> (-1)*sum(V[j,k] * Basefunc[j](x) for j in 1:n)
end

#=====================================================#
#--------------- Loop Principal------------------------- 

function despejar_frecuencia(L::Vector,A::Vector, a::Number, m::Number,tono::Int) #L:=autoval R A:= autoval S
    """Encuentra el  valor de w  para un  autovalor del  operador S"""
    n = length(A)
    soluciones = Vector{ComplexF64}(undef, n)
    for  i in 1:n
        p = Polynomial([(A[i] - L[i]), (-2*a*m), (a^2)])   # estes  es E(w;A)
        w = roots(p)                                       # aqui hago que tienda a cero
        soluciones[i] = w[tono]
    end
    return soluciones
end

function calculo_error_abs(v_0::Number,v_1::Number)
    n_re = abs(real(v_1) - real(v_0))
    n_im = abs(imag(v_1) - imag(v_0))
    return sqrt(n_re^2 + n_im^2)
end

function actualizar_convergencia(iteraciones::Vector,error_A_array::Vector,error_w_array::Vector,N::Number,error_w::Number,error_A::Number)
    iteraciones[N] = N
    error_w_array[N] = error_w
    error_A_array[N] = error_A
end


function iteracion_newton(a::Number,m::Int,M::Number,n::Int,w_inicial::ComplexF64,
                             s::Int,N_max::Int,error_max::Float64,tono::Int,l::Int)
    iteraciones = zeros(Int, N_max)
    error_w_array = zeros(Float64, N_max)
    error_A_array = zeros(Float64, N_max)
    error_w = 0.0
    error_A = 0.0
    N=1; error=1e2;A_inicial=0

    while N<=N_max && error >=error_max 

        println("\n iteración N=$N \n")
        A_ang,V_ang=sol_espectral(operador_angular, w_inicial,:ii, a, s, m , M, n) #A:= auto  valores  V:=auto vectores
        println("A_ ang= $A_ang \n")
       
        U_rad,V_rad=sol_espectral(operador_radial, w_inicial, :i, a, s, m , M, n) 
        println("U_RAD N=$U_rad \n")
       
        w_act= despejar_frecuencia(U_rad,A_ang, a, m,tono)
        println("frecuencias =$w_act \n")


        error_w=calculo_error_abs(w_inicial,w_act[l])
        error_A=calculo_error_abs(A_inicial,A_ang[l])

        actualizar_convergencia(iteraciones,error_A_array,error_w_array,N,error_w,error_A)

        w_inicial=w_act[l]
        A_inicial=A_ang[l]
        N+=1 ; error=error_w + error_A
        

    end

    if N > N_max
        println("Máximo número de iteraciones alcanzado")
    end
    if error < error_max
        println("Convergencia alcanzada con error en la frecuencia=$error_w  y el error en A=$error_A")
    end
    
    return  iteraciones,error_A_array,error_w_array, w_inicial,A_inicial, V_ang,V_rad
end

#=====================================================#
#-------------------- Graficas------------------------- 


function PL_m1(x, l)  # comparación  teorica  con la parte  angular 
    """Polinomios de Legendre asociados P_l^1(x)  (m=1)"""
    
    if l == 0 # Derivada analítica de cada P_l(x)
        return 0.0
    elseif l == 1
        dP = 1.0
    elseif l == 2
        dP = 3x
    elseif l == 3
        dP = 0.5 * (15x^2 - 3)
    elseif l == 4
        dP = (1/8) * (140x^3 - 60x)
    elseif l == 5
        dP = (1/8) * (315x^4 - 210x^2 + 15)
    else
        error("Solo definido hasta l=5 por ahora")
    end
    return (-1)*sqrt(1 - x^2) * dP     # P_l^1(x)
end

function R_ansatz(H::Function, M::Number, a::Number, w::Number, m::Number)  # comparación  teorica  con la parte  radial 
    # los horizontes
    rmas= M + sqrt(M^2-a^2)
    rmenos= M - sqrt(M^2-a^2)
    R_pm= rmas-rmenos 
    i = 0 + 1im 
    
    #elemntos del cambio de  variable z=(r-r+)/(r-r-)
    r= z-> (z*rmenos-rmas)/(z-1)
    
    #elementos del ansatz ecu.5.20
    sigmaMas=(2*w*M*rmas-m*a)/(R_pm)
    sigmaMenos=(2*w*M*rmenos-m*a)/(R_pm)
    
    # se toma el  signo positivo en el pm
    zeta=i*w
    xi=i*sigmaMas  
    eta=-i*sigmaMenos
    # se toma el signo negativo en el pm 
    #zeta=-i*w
    #xi=-s-i*sigmaMas  
    #eta=-s+i*sigmaMenos

    return z-> (r(z)-rmas)^xi * (r(z)-rmenos)^eta * exp(zeta*r(z)) * H(z)

end



function grafica5_modos_comparativa(V_exp::Matrix,U_exp::Vector,PL_teor::Function,tipo::Symbol, parte::String,  M::Number, a::Number, w::Number, m::Number)
   
    
    n_modos = 5
    p = plot(layout = (n_modos, 1), size=(700, 1000))

    for k in 1:n_modos
        Yk = Y_k_aproximado(V_exp, tipo, k)
        
        if parte=="ang"
            x = range(-1, 1, length=400)

            plot!(p[k], x, [real(Yk(xi)) for xi in x],
            label = latexstring("A_{1\\ell} = ", round(real(U_exp[k]), digits=2)),
            lw = 2)

            plot!(p[k], x, [PL_teor(xi, k) for xi in x],
                ls = :dash, lw = 2, color = :black,
                label = L"P_{\ell}^{1}(x) ")

            xlabel!(p[k], L"x = \cos(\theta)")
            ylabel!(p[k], L"S_{\ell 1}(x)")

        elseif  parte== "rad"
            x = range(0, 0.95, length=400)

            plot!(p[k], x, [real(Yk(xi)) for xi in x],
            label = latexstring("H_{1\\ell}(z)"),
            lw = 2)

            R=R_ansatz(Yk, M, a, w, m)
            plot!(p[k], x, [real(R(xi))  for xi in x],
                ls = :dash, lw = 2, color = :black,
                label = L"R_{\ell}^{1}(z)")
            xlabel!(p[k], L"z =\frac{r-r_{+}}{r-r_{-}}")
            ylabel!(p[k], L".")
        end

        
        title!(p[k], "Modo ℓ = $k")
    end

    plot!(p, legend = :bottomleft)
    display(p)

end

function grafica_convergencia(iteraciones::Vector,error_A_array::Vector,error_w_array::Vector, N::Int)

    rango = 1:N        # tomar solo las iteraciones válidas

    plot(iteraciones[rango], error_w_array[rango],
         lw = 2, label = L"\epsilon_{\omega}",
         xlabel = "Iteración", ylabel = "Error", yscale = :log10)

    plot!(iteraciones[rango], error_A_array[rango],
          lw = 2, ls = :dash, label = L"\epsilon_{A}",
          legend = :bottomleft)

    title!("Convergencia del método")
    plot!(legend = :topright)
end




### Fín del  módulo 
end
