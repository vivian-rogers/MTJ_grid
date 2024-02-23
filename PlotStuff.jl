using Plots
using SparseArrays
using LinearAlgebra
using Plots
import Plots
using IndirectArrays, FileIO, Colors

function Mvec_to_matrices(p::NamedTuple, M::Matrix{Float64},dim::Int=3)
    nT = size(M)[2]
    Mzgrid = zeros((p.nx,p.ny,nT))
    display(size(M))
    for it = 1:nT
        Mzs_at_t = M[dim:3:(p.nx*p.ny*3),it]
        #display(sizeof(Mzs_at_t))
        #Mzs_at_t = M[(3):3:(p.nx*p.ny),it]
        Mzgrid[:,:,it] = reshape(Mzs_at_t,(p.nx,p.ny))
    end
    return Mzgrid
end

#commented because @animate for some reason not recognized
#=function animate_Mzvals(p::NamedTuple, Mzgrid::Array{Float64, 3})
    nT = size(Mzgrid)[3]
    anim = @animate for it ∈ 1:nT
        heatmap(Mzgrid[:,:,it]',clim=(-1,1), cmap=:coolwarm, title="$it")
    end
    return gif(anim, "anim_fps15.gif", fps = 30)
end=#

function sum_Mz_time(p::NamedTuple, Mzgrid::Array{Float64,3})
    nT = size(Mzgrid)[3]
    avgMzs = [mean(Mzgrid[:,:,it]) for it = 1:nT]
    stdMzs = [std(Mzgrid[:,:,it]) for it = 1:nT]
    plot(avgMzs, label="Avg M");
    plot!(stdMzs, ylim=(-1.1,1.1), label="Std(M)")
end

function final_Mzvals(p::NamedTuple, Mzgrid::Array{Float64, 3})
    nT = size(Mzgrid)[3]
    return heatmap(Mzgrid[:,:,nT],clim=(-1,1), cmap=:coolwarm, title="finaltime")
end

function isolate_device(p::NamedTuple, M::Array, ivec::Vector{Int})
    nT = size(M)[2]
    M_over_t = zeros(3,nT)
    idevice = ivec_to_index(p, ivec)
    for it = 1:nT
        M_over_t[:,it] = M[(3*idevice-2):(3*idevice),it]
    end
    return M_over_t
end

function plot_single_device(p::NamedTuple, M::Array, ivec::Vector{Int})
    nT = size(M)[2]
    M_t = isolate_device(p,M,ivec)
    plot(M_t[1,:], M_t[2,:], M_t[3,:], xlim=(-1.1,1.1), ylim=(-1.1,1.1), zlim=(-1.1,1.1), ms=5)
    #anim = @animate for it ∈ 1:nT
    #    scatter(M_t[1,it], M_t[2,it], M_t[3,it], title="$it",ms=5)
    #end
   # return gif(anim, "anim_singledevice_fps15.gif", fps = 120)
end


function fullColorAnimation(p::NamedTuple, Mvals::Matrix{Float64}, name::String="animation.gif")
    nT = size(Mvals)[2]
    Mgrid = zeros((3,p.nx,p.ny,nT))
    Mgrid[1,:,:,:] = (Mvec_to_matrices(p,Mvals,1))
    Mgrid[2,:,:,:] = (Mvec_to_matrices(p,Mvals,2))
    Mgrid[3,:,:,:] = (Mvec_to_matrices(p,Mvals,3).+1)/2
    H = (atan.(Mgrid[1,:,:,:], Mgrid[2,:,:,:]).+π)*(180/π)
    S = ones(size(Mgrid[1,:,:,:]))
    L = Mgrid[3,:,:,:]
    #colorarray = IndirectArray(RGB.(Mgrid), [colorant"red", colorant"green", colorant"blue"])
    save(name, HSL.(H,S,L); fps=40)
end

function plotSpectrum(p::NamedTuple, Mvals::Matrix{Float64}, name::String="spectrum.png", dim::Int=3)
    Mzgrid = Mvec_to_matrices(p,Mvals,dim)
    nT = size(Mzgrid)[3]
    Mzvals = Mzgrid[:,:,nT]
    eig = eigvals(sign.(Mzvals))
    fig = plot(heatmap(Mzgrid[:,:,nT],clim=(-1,1), cmap=:coolwarm, title="finaltime"), 
    scatter(real.(eig), imag.(eig), xlabel="Re(λ)", ylabel="Im(λ)"), size=(900,450))
    savefig(fig, name)
end