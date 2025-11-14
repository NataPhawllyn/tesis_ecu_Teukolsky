# Comprobamos si el paquete ya está instalado antes de agregarlo
function add_package_if_missing(pkg)
    if !(pkg in keys(Pkg.installed()))
        println("Instalando el paquete: $pkg")
        Pkg.add(pkg)
    else
        println("El paquete $pkg ya está instalado.")
    end
end

# Lista de paquetes a instalar
packages = [
    "LinearAlgebra", "SparseArrays", "NLsolve", "ForwardDiff",
    "DifferentialEquations", "Plots", "PlotlyJS", "LaTeXStrings",
    "GLMakie", "CairoMakie", "DelimitedFiles", "CSV", "DataFrames", "Printf"
]

# Instalar los paquetes si no están presentes
for pkg in packages
    add_package_if_missing(pkg)
end

# Eliminar paquetes innecesarios
Pkg.rm("GIFs")
Pkg.rm("PyPlot")
Pkg.rm(["PyCall", "SymPy", "Makie"])
Pkg.gc()

# Paquetes adicionales
add_package_if_missing("QuadGK")
add_package_if_missing("SpecialPolynomials")
add_package_if_missing("Polynomials")
add_package_if_missing("Roots")

# Precompilación de paquetes
Pkg.precompile()

println("Instalación completada.")




# import Pkg
# Pkg.add([
#     "LinearAlgebra",
#     "SparseArrays",
#     "NLsolve",
#     "ForwardDiff",
#     "DifferentialEquations",
#     "Plots",
#     "PlotlyJS",
#     "LaTeXStrings",
#     "GLMakie",
#     "CairoMakie",
#     "DelimitedFiles",
#     "CSV",
#     "DataFrames",
#     "Printf"
# ])
# Pkg.rm("GIFs")
# Pkg.rm("PyPlot")
# Pkg.rm(["PyCall", "SymPy", "Makie"])
# Pkg.gc()
# Pkg.add("QuadGK")
# Pkg.add("GLMakie")
# Pkg.precompile()
# Pkg.add("SpecialPolynomials")
# Pkg.add("Polynomials")
# using Pkg; Pkg.add("Roots")
# using Pkg; Pkg.add("NLsolve")