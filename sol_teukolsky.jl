using ForwardDiff, QuadGK
using LinearAlgebra
using SpecialPolynomials
using Polynomials
using Plots, LaTeXStrings
using .Utilidades 


#================================================#
#--------- Ejemplo  limite de Schwarzschild-------
# a=0 ; w = 4+12im ; s = 0; m = 1 ; M=1; n=5 #tamaño de la  base 

#...Parte Angular
# tipo=:ii  # Polinomios TIPO U
# U_ang,V_ang=sol_espectral(operador_angular, w, tipo, a, s, m , M, n)
# P=[2,6,12,20,30]  # autovalores teoricos 
# println("autovalores parte angular :\n",U_ang)
# grafica5_modos_comparativa(V_ang,U_ang,PL_m1,tipo,"ang",M, a, w, m) # encontrado  vs polinomios de Laguerre 

# ...Parte Radial 
# tipo=:i  #Polinomios TIPO T
# U_rad,V_rad=sol_espectral(operador_radial, w, tipo, a, s, m , M, n)
# println("autovalores parte  radial:\n",(U_rad))
# grafica5_modos_comparativa(V_rad,U_rad,R_ansatz,tipo,"rad", M, a, w, m) #  H  vs R  encontrados


#-------------Ejemplo Bertti - Metodo Leaver--------#
# L=3; a=0 ; s = 2; m = 1 ; M=1;
# Alm=L*(L+1) ; w = 0.005+0.02im    # valor esperado  0.59944-0.09270im  Mw, 1.19889+0.185406im 2Mw 
# w_nuevo= despejar_omega(w, Alm, a, s, m , M,"Shwa")
# println(w_nuevo)  # unidades Mw se divide en 2

#================================================#
#============          KERR        ==============#

#================================================#
#------------------ LOOP GENERAL ------------------

a=0 ; m = 1 ; M=1; n=5   # constantes fijas,  n:= tamaño de la  base 
w_inicial = 0.6-2im ; s = 0  #+2, -2  para ondas  gravitacionales
N_max= Int(1e2); error_max=1e-12  #Numero de iteraciones máximas y  termino de  tolerancia
tono=1; l=3  #se elije el tono de w y la linea de modos |m|<=l

iteraciones,error_A_array,error_w_array, w_final,A_final, V_ang=iteracion_newton(a,m,M,n,w_inicial,s,N_max,error_max,l)

println("frecuencia encontrada= $(w_final/2), A  encontrado $A_final")
# grafica_convergencia(iteraciones,error_A_array,error_w_array,N_max)
