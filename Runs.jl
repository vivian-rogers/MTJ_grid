include("./Constants.jl")
include("./PlotStuff.jl")
include("./Solver.jl")


# all values in Kg-m-s units unless otherwise specified

parameters = (Ms = 1.5E6, Δz = 1.5*nm, Δx = 20*nm, Δy = 20*nm, μᵣ = 0.5, Δt = 2E-12, t_final = 5E-9, T=10,
nx = 15, ny = 15, ax = 30*nm, ay=30*nm, α = 0.1, Ku = 1.0E6, C=1, B_app = (B_app(t::Float64)= 0.01*ẑ));

Mvals, tvals = heun_integrate(parameters)
fullColorAnimation(parameters, Mvals, "test.gif")