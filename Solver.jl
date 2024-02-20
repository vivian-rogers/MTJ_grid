using LinearAlgebra
using SparseArrays
using Random, Distributions
using Distributed
using Statistics
using BlockDiagonals

# some useful definitions here
Random.seed!(123)
g = Normal();
⊗(A,B) = kron(A,B)


# for the cross product A×B, expresses (A×) as a matrix 
crossmat(A::Vector{Float64}) = [0 -A[3] A[2]; A[3] 0 -A[1]; -A[2] A[1] 0]

# gets the devices' device scalar index as a function of the [ix; iy] index vector
function ivec_to_index(p::NamedTuple, ivec::Vector{Int})
    ix = ivec[1]; iy = ivec[2]
    return ix + 1  + iy*p.nx
end

# gets the position from a device's index vector
function ivec_to_position(p::NamedTuple, ivec::Vector{Int})
    return ivec[1]*p.ax*x̂ + ivec[2]*p.ay*ŷ
end

# gets the magnetization from a devices' device index
function idevice_to_M_device(M::Vector{Float64},idevice::Int)
    return M[(3*idevice-2):(3*idevice)]
end

# gets the magnetization from a devices' index vector
function ivec_to_M_device(p::NamedTuple, M::Vector{Float64},ivec::Vector{Int})
    idevice = ivec_to_index(p, ivec)
    return idevice_to_M_device(M, idevice)
end

# calculate the stray fields emanating from one device
function H_from_deviceA_on_deviceB(p::NamedTuple, ivecA::Vector{Int}, ivecB::Vector{Int}, M::Vector{Float64})
    rA = ivec_to_position(p, ivecA)
    rB = ivec_to_position(p, ivecB)
    ΔR = rB-rA; R = norm(ΔR); r̂ = ΔR/R
    m̂A = ivec_to_M_device(p, M, ivecA)
    H = p.m*(3*r̂*(m̂A⋅r̂) - m̂A)/(4*π*R^3)
    return H
end

# this will generate the matrix that is used in the statement of dψ/dt = A(t)ψ for the solver
function Agen(p::NamedTuple, M::Vector{Float64}, t::Float64)
    blockMatrices = Vector{Matrix{Float64}}(undef,p.nx*p.ny)
    # you could put a @Threads.threads macro here
    for ix = 0:(p.nx-1)
        for iy = 0:(p.ny-1)
            ivec = [ix;iy]
            idevice = ivec_to_index(p,ivec)
            m̂_device = idevice_to_M_device(M,idevice)
            # now let's build the Heff 
            random_vector = rand(g,3); random_vector = random_vector/norm(random_vector);
            Htherm = random_vector*√(2*μ₀*p.α*p.T*kB/(p.Ms*μ₀*p.Δt*p.ΔV*abs(γ)))/μ₀; 
            H_anisotropy = (2*p.Ku/(μ₀*p.Ms))*m̂_device[3]*ẑ
            
            # now throw it all together
            Heff = Htherm + H_anisotropy + p.B_app(t)/μ₀

            # can uncommend this if you wish to implement ferromagnetic-like exchange interactions
            #=Δivec = zeros(2)
            for Δivec ∈ [[-1,0],[1,0],[0,-1],[0,1]]
                ivecB = ivec + Δivec
                if(ivecB[1]<=(p.nx-1) && ivecB[1]>=0 && ivecB[2]<=(p.ny-1) && ivecB[2]>=0)
                    Heff += p.C*H_from_deviceA_on_deviceB(p, ivecB, ivec, M)
                end
            end=#
            for Δix ∈ [-1,1]
                for Δiy ∈ [-1,1]
                    ivecB = ivec + [Δix; Δiy]
                    if(ivecB[1]<=(p.nx-1) && ivecB[1]>=0 && ivecB[2]<=(p.ny-1) && ivecB[2]>=0)
                        Heff += p.C*H_from_deviceA_on_deviceB(p, ivecB, ivec, M)
                    end
                end
            end
            #println("ivec = $ivec, Heff = $Heff")
            # THIS IS THE LL EQUATION 
            A_device = γ*crossmat(μ₀*Heff)  + p.λ*crossmat(cross(m̂_device,μ₀*Heff))
            # and we throw it onto the big diagonal matrix 
            blockMatrices[idevice] = copy(A_device)
            #=for icol = 1:3
                for irow = 1:3
                    push!(colvals, 3*idevice-3+icol)
                    push!(rowvals, 3*idevice-3+irow)
                    push!(zvals, A_device[irow, icol])
                end
            end=#
        end
    end
    return BlockDiagonal(blockMatrices)
    #return sparse(rowvals, colvals, zvals)
    #return sparse(rowvals, colvals, zvals)
end

# calculates more parameters as a function of your input parameters
function more_parameters(p::NamedTuple)
    λ = p.α*γ/(μ₀*p.Ms)
    ΔV = p.Δy*p.Δx*p.Δz
    m = p.Ms*ΔV
    tvals = collect([t for t = 0:p.Δt:p.t_final]); nT = size(tvals)[1]
    newparams = merge(p, (λ = λ, ΔV = ΔV, m = m, nT = nT))
    display(newparams)
    return newparams
end

# performs the Heun integration method to simulate your system for N timesteps
function heun_integrate(p::NamedTuple)
    p = more_parameters(p)
    tvals = collect([t for t = 0:p.Δt:p.t_final]); nT = size(tvals)[1]
    println("Integrating over $nT timesteps")
    nD = p.nx*p.ny
    Mvals = zeros(nD*3,nT)
    # initialize it with something 
    M₀ = ones(nD)⊗[1;0;0]
    #=M₀ = rand(g, 3*nD)
    # normalize each m\hat
    for iD ∈ 1:nD
        m̂ = M₀[(1+3*iD-3):(3*iD)]
        M₀[(1+3*iD-3):(3*iD)] = m̂/norm(m̂)
    end=#
    for it ∈ eachindex(tvals)
        #println("t = $(tvals[it])")
        t = tvals[it]
        f₀ = Agen(p,M₀,t)*M₀
        M̃₁ = M₀ + p.Δt*f₀
        for iD ∈ 1:nD
            m̂ = M̃₁[(1+3*iD-3):(3*iD)]
            M̃₁[(1+3*iD-3):(3*iD)] = m̂/norm(m̂)
        end
        #M̃₁ = M₀ + Δt*Agen(p,M₀)*M₀
        M₁ = M₀ + (p.Δt/2)*(f₀+Agen(p,M̃₁,t)*M̃₁)
        # and normalize
        for iD ∈ 1:nD
            m̂ = M₁[(1+3*iD-3):(3*iD)]
            M₁[(1+3*iD-3):(3*iD)] = m̂/norm(m̂)
        end
        Mvals[:,it] = M₁
        M₀ = deepcopy(M₁)
    end
    println("Done! Returning magnetization values")
    return Mvals, tvals
end
