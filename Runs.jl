include("./Constants.jl")
include("./PlotStuff.jl")
include("./Solver.jl")


# all values in Kg-m-s units unless otherwise specified

parameters = (Ms = 1.5E6, Δz = 1.5*nm, Δx = 50*nm, Δy = 50*nm, μᵣ = 0.5, Δt = 2.0E-12, t_final = 25.0E-9, T=5000, init="x",
nx = 30, ny = 30, ax = 60*nm, ay=60*nm, α = 0.01, Ku = 0.10E6*ẑ, C=1, B_app = (B_app(t::Float64)= 0.000*ẑ));

@time Mvals, tvals = heun_integrate(parameters)
plotSpectrum(parameters, Mvals, "30x30spectrum.png",3)
fullColorAnimation(parameters, Mvals, "30x30test.gif")