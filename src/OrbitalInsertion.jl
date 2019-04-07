__precompile__()

module OrbitalInsertion

using DifferentialEquations

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

export search, createConstants, createPar

include("funcs_insert.jl")

end # module
