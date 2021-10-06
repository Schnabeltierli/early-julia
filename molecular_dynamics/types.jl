mutable struct Particle
    pos::Vector{Float64}
    vel::Vector{Float64}
    acc::Vector{Float64}
    μ::Float64
    σ::Float64
    ϵ::Float64
end
