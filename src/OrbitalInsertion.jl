__precompile__()

module OrbitalInsertion

using DifferentialEquations, ForwardDiff, DiffResults, Printf, LinearAlgebra, NLopt, Roots, Interpolations, Plots

abstract type System end

struct EMsystem <:System end
struct SEsystem <:System end

struct ProblemInfo{T<:Float64}
    mu  ::T
    rp  ::T
    rs  ::T
    r12 ::T
    w   ::T
    mup ::T
    x1  ::T
    x2  ::T
end
struct ProblemParameters{T<:Float64,J<:Int64,K<:Bool}
    tspan   ::T
    hgoal   ::T
    hmmin   ::T
    hmmax   ::T
    sType   ::K
    points  ::J
    stop    ::K
    maxt0   ::T
    rest    ::Tuple{J,J}
    lim     ::T
    cTol    ::T
    tRest   ::K
    odeAlg  ::K
    atol    ::T
    rtol    ::T
    retorno ::K
    detMin  ::T
    x0      ::Array{T,1}
    y       ::Array{T,1}
    optTol  ::T
    neval   ::J
    col     ::K
end
mutable struct SystemInfo{T<:Float64}
    mu  ::T
    Rp  ::T
    Rs  ::T
    r12 ::T
    m1  ::T
    G   ::T
end

export createConstants, createSystem, createPar, search, filterGlobalSol, point, plotSols
export EMsystem, SEsystem
export ProblemInfo, ProblemParameters, SystemInfo, ProblemParameters

include("funcs_insert.jl")
include("systems.jl")
include("restrictions.jl")

end # module
