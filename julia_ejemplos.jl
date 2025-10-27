using Printf
using LinearAlgebra
using Plots

sumsq(x,y)=x^2+y^2;
N=12
x=0
for  i =1:N
    if  sumsq(rand(),rand())<1.0
        global x+=1
    end
end

@printf "Hola  mundo  %d es %8.5f\n " N*2 4.0*(x/N)

A=[1 -2 2; 1 -1 2; -1 1 1]
U=eigvals(A)
V=eigvecs(A)

@printf("Los autovalores son:\n%s\n\n", U)
@printf("Los autovectores son:\n%s\n", V)

yy = 2.0*randn(100,3)

plot(yy, layout=grid(3,1, heights=[0.4,0.3,0.3]))


