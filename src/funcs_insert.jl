"""
Description createSystem
"""
function createSystem(sys::System,G=6.673e-11::Float64)
    if      typeof(sys) == EMsystem
        return SystemInfo(EM.mu,EM.Rp,EM.Rs,EM.r12,EM.m1,Float64(G))
    elseif  typeof(sys) == SEsystem
        return SystemInfo(SE.mu,SE.Rp,SE.Rs,SE.r12,SE.m1,Float64(G))
    end
end
"""
Description createConstants
"""
function createConstants(sysData::SystemInfo)
    G       = sysData.G # Gravity constant      [m3/Kg/s2]
    ### Make sure that the arguments are Floats
    mu      = sysData.mu
    if mu > 1 || mu < 0
        error("mu parameter has to be 0<mu<1")
    end
    Rp      = sysData.Rp
    Rs      = sysData.Rs
    r12     = sysData.r12
    m1      = sysData.m1
    ###
    m2      = m1*mu/(1-mu)
    # Rotation of reference system [rad/s]
    w       =  sqrt(G*(m1+m2)/r12^3)
    # Gravity constant of Earth [m^3/s^2]
    mup     = G*m1
    return ProblemInfo(mu,Rp/r12,Rs/r12,r12,w,mup/r12^3/w^2,-mu,1-mu)
end
"""
Description createPar
"""
function createPar(prob_info;   tspan = -2.0, hgoal = 300, hmmin = 100 ,
                        hmmax = 20e3, sType = "global", points = 20, x0 = [],
                        maxt0 = -1.0, rest = "basic", lim = 2000,
                        cTol = 1e-4, tRest = "ineq", alg = "lowOrder",
                        atol = 1e-7, rtol = 1e-6, retorno = true, detMin = true, y = [], optTol = 1e-8, neval = 400)

    if length(x0) == 0
        println("There has to be an initial guess ('x0' argument)")
        return nothing
    elseif length(x0) == 3
        stop = true
    elseif length(x0) == 4 || length(x0) == 6
        detMin  = 1.0
        stop    = true
    elseif length(x0) == 5 || length(x0) == 7
        stop = false
    end

    if length(x0) == 4 || length(x0) == 5
        col = true
    else
        col = false
    end

    if sType == "global"
        type = true
    elseif sType == "local"
        type = false
    else
        println("The search type has to be either 'global' or 'local' ('sType' argument)")
        return nothing
    end

    if tRest == "ineq"
        restType = true
    elseif tRest == "mixed"
        restType = false
    else
        println("The restriction type has to be 'ineq' or 'mixed' ('tRest' argument)")
        return nothing
    end
    restriction = dictRest[rest]

    if alg == "lowOrder"
        odeAlg = true
    elseif alg == "highOrder"
        odeAlg = false
    else
        println("Wrong ODE algorithm, choose between 'lowOrder' or 'highOrder'")
        return nothing
    end

    if y == []
        println("There has to be an insertion point ('y' argument)")
        return nothing
    elseif length(y) > 6 || length(y) < 6
        println("The insertion point ('y' argument) has to be of length 6")
        return nothing
    end

    println("--- Problem parameters created correctly ---")
    return ProblemParameters(Float64(tspan),Float64(hgoal*1e3/prob_info.r12),Float64(hmmin*1e3/prob_info.r12),Float64(hmmax*1e3/prob_info.r12),type,Int64(points),stop,Float64(maxt0),
                            restriction,Float64(lim/prob_info.r12/prob_info.w),Float64(cTol),restType,odeAlg,Float64(atol),Float64(rtol),
                            retorno,Float64(detMin),convert.(Float64,x0),convert.(Float64,y),optTol,Int64(neval),col)
end
"""
Description dineq!
"""
function dineq!(dy,y,p,t,mu::Float64)
    #mu  = p[1]
    xp, yp, zp, vx, vy, vz = y
    ro1 = sqrt((xp+mu)^2      +yp^2 +zp^2)
    ro2 = sqrt((xp-(1-mu))^2  +yp^2 +zp^2)
    dy[1] = vx
    dy[2] = vy
    dy[3] = vz
    dy[4] = 2*vy+xp-(1-mu)*(xp+mu)/ro1^3-mu*(xp-(1-mu))/ro2^3
    dy[5] = -2*vx+yp*(1 - (1-mu)/ro1^3-mu/ro2^3)
    dy[6] = -zp*((1-mu)/ro1^3+mu/ro2^3)
end
"""
Description calc_data
"""
function calc_data(x::Array{Float64,1},par::ProblemParameters, prob_info::ProblemInfo, PROB::Union{Tuple{ODEProblem,ODEProblem},Tuple{ODEProblem}})
    x1      = prob_info.x1      ::Float64
    x2      = prob_info.x2      ::Float64
    rearth  = prob_info.rp      ::Float64
    rmoon   = prob_info.rs      ::Float64
    prob    = PROB[1]
    y       = prob.u0           ::Array{Float64,1}

    ny0 = y - vcat(zeros(3),x[1:3])

    if par.stop
        if length(x) == 3
            intProb          = remake(prob, u0 = ny0 )
        elseif length(x) == 4
            intProb          = remake(prob, u0 = ny0 , p = [x[4], 0.0, 0.0, 0.0, par.detMin])
        else
            intProb          = remake(prob, u0 = ny0 , p = [x[4], x[5], x[6], 0.0, par.detMin])
        end
        tstop = [0.0]
    else
        if length(x) == 5
            intProb          = remake(prob, u0 = ny0 , p = [x[4], 0.0, 0.0, x[5], par.detMin])
        elseif length(x) == 7
            intProb          = remake(prob, u0 = ny0 , p = [x[4], x[5], x[6], x[7], par.detMin])
        end
        tstop = [x[end]]
    end

    if par.odeAlg
        sol     = solve(intProb, DP5(), abstol = par.atol, reltol= par.rtol, tstops= tstop, dense=false, save_everystep=false)
    else
        sol     = solve(intProb, Vern7(), abstol = par.atol, reltol= par.rtol, tstops= tstop, dense=false, save_everystep=false)
    end

    ts, tp, ps, pp, ImpF, ImpS, ImpI, tImpF, tImpS, tImpI = calcOuts(x,sol,par,prob_info)

    if par.retorno
        probR = PROB[2]
        if length(x) == 3
            Y0 = ImpF
        elseif length(sol.u) > 2
            Y0 = ImpS
        else
            Y0 = ImpF
        end
        prob2   = remake(probR, u0 = Y0 , tspan = (0,2.5) )
        if par.odeAlg
            solr    = solve(prob2, DP5(), abstol = par.atol, reltol= par.rtol, dense=false, save_everystep=false)
        else
            solr    = solve(prob2, Vern7(), abstol = par.atol, reltol= par.rtol, dense=false, save_everystep=false)
        end
        eRet       = calcDists(solr.u[end], x1,rearth)
    else
        eRet       = 0.0
    end


    ds, dp   = calcDists(ps, x2, rmoon), calcDists(pp, x1, rearth) # distances to primaries in the minimum distance points

    fimp = fImp(pp,prob_info)

    return ds, dp, ts, tp, eRet, ImpF, ImpS, ImpI, tImpF, tImpS, tImpI, fimp
end
# Add method to calc_data
function calc_data(x::Vector{T},par::ProblemParameters, prob_info::ProblemInfo, PROB::Tuple{ODEProblem}) where {T<:ForwardDiff.Dual}
    x1      = prob_info.x1      ::Float64
    x2      = prob_info.x2      ::Float64
    rearth  = prob_info.rp      ::Float64
    rmoon   = prob_info.rs      ::Float64
    prob    = PROB[1]

    xType = eltype(x)
    if par.stop
        if length(x) == 3
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x), tspan = convert.(xType, prob.tspan) )
        elseif length(x) == 4
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], 0.0, 0.0, 0.0, par.detMin]) )
        else
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], x[5], x[6], 0.0, par.detMin]) )
        end
        tstop = convert.(xType,[0.0])
    else
        if length(x) == 5
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], 0.0, 0.0, x[5], par.detMin]) )
        elseif length(x) == 7
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], x[5], x[6], x[7], par.detMin]) )
        end
        tstop = x[end]
    end

    if par.odeAlg
        sol     = solve(intProb, DP5(), abstol = par.atol, reltol= par.rtol, tstops=tstop, dense=false, save_everystep=false)
    else
        sol     = solve(intProb, Vern7(), abstol = par.atol, reltol= par.rtol, tstops=tstop, dense=false, save_everystep=false)
    end

    ts, tp, ps, pp = calcOuts(x,sol,par,prob_info,PROB)

    #ds, dp   = calcDists(ps, x2, rmoon), calcDists(pp, x1, rearth) # distances to primaries in the minimum distance points
    ds = calcDists(ps, x2, rmoon)
    dp = calcDists(pp, x1, rearth)

    return ds, dp, ts, tp
end
function calc_data(x::Vector{T},par::ProblemParameters, prob_info::ProblemInfo, PROB::Tuple{ODEProblem,ODEProblem}) where {T<:ForwardDiff.Dual}
    x1      = prob_info.x1      ::Float64
    x2      = prob_info.x2      ::Float64
    rearth  = prob_info.rp      ::Float64
    rmoon   = prob_info.rs      ::Float64
    prob    = PROB[1]
    probR   = PROB[2]

    xType = eltype(x)
    if par.stop
        if length(x) == 3
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x), tspan = convert.(xType, prob.tspan) )
        elseif length(x) == 4
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], 0.0, 0.0, 0.0, par.detMin]) )
        else
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], x[5], x[6], 0.0, par.detMin]) )
        end
        tstop = convert.(xType,[0.0])
    else
        if length(x) == 5
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], 0.0, 0.0, x[5], par.detMin]) )
        elseif length(x) == 7
            intProb          = remake(prob, u0 = prob.u0 - vcat(zeros(3),x[1:3]), tspan = convert.(xType, prob.tspan), p = convert.(xType,[x[4], x[5], x[6], x[7], par.detMin]) )
        end
        tstop = x[end]
    end

    if par.odeAlg
        sol     = solve(intProb, DP5(), abstol = par.atol, reltol= par.rtol, tstops=tstop, dense=false, save_everystep=false)
    else
        sol     = solve(intProb, Vern7(), abstol = par.atol, reltol= par.rtol, tstops=tstop, dense=false, save_everystep=false)
    end

    ts, tp, ps, pp, IMP = calcOuts(x,sol,par,prob_info,PROB)

    prob2   = remake(probR, u0 = IMP , tspan = convert.(xType,(0.0,2.5)) )
    if par.odeAlg
        solr    = solve(prob2, DP5(), abstol = par.atol, reltol= par.rtol, dense=false, save_everystep=false, callback = callRet )
    else
        solr    = solve(prob2, Vern7(), abstol = par.atol, reltol= par.rtol, dense=false, save_everystep=false, callback = callRet )
    end
    eRet       = calcDists(solr.u[end], x1,rearth)

    #ds, dp   = calcDists(ps, x2, rmoon), calcDists(pp, x1, rearth) # distances to primaries in the minimum distance points
    ds = calcDists(ps, x2, rmoon)
    dp = calcDists(pp, x1, rearth)

    return ds, dp, ts, tp, eRet
end
"""
Description calcOuts
"""
function calcOuts(x::Array{Float64,1},sol::ODESolution,par::ProblemParameters,prob_info::ProblemInfo)
    ts = sol.t[1] # Instant of close approach to secondary body
    tp = sol.t[end] # Instant of close approach to primary body
    ps = sol.u[1] # State of close approach to secondary body
    pp = sol.u[end] # State of close approach to primary body

    ImpF    = sol.u[1]
    tImpF   = sol.t[1]
    ImpI    = sol.u[end]
    tImpI   = sol.t[end]

    ImpS    = sol.u[2]
    tImpS   = sol.t[2]
    if length(x) == 3
        if length(sol.t) == 3
            ts = sol.t[2]
            ps = sol.u[2]
        end
    else
        if par.stop
            if length(sol.t) == 3
                ts = sol.t[2]
                ps = sol.u[2]
            end
        else
            if length(sol.t) >= 4
                tImp        = x[end]

                indexImp = findfirst( x -> sol.t[x] == tImp, sol.t)
                ImpS  = sol.u[indexImp]
                tImpS   = sol.t[indexImp]

                Mclose = [ (calcDists(sol.u[i],prob_info.x2,prob_info.rs),i) for i in 2:(length(sol.t)-1) if i != indexImp ]
                sort!(Mclose,by = x -> x[1])
                ts = sol.t[Mclose[1][2]]
                ps = sol.u[Mclose[1][2]]
            end
        end
    end
    return ts, tp, ps, pp, ImpF, ImpS, ImpI, tImpF, tImpS, tImpI
end
function calcOuts(x::Vector{T},sol::ODESolution,par::ProblemParameters,prob_info::ProblemInfo,PROB::Tuple{ODEProblem}) where {T<:ForwardDiff.Dual}
    ts = sol.t[1] # Instant of close approach to secondary body
    tp = sol.t[end] # Instant of close approach to primary body
    ps = sol.u[1] # State of close approach to secondary body
    pp = sol.u[end] # State of close approach to primary body

    if length(x) == 3
        if length(sol.t) == 3
            ts = sol.t[2]
            ps = sol.u[2]
        end
    else
        if par.stop
            if length(sol.t) == 3
                ts = sol.t[2]
                ps = sol.u[2]
            end
        else
            if length(sol.t) >= 4
                tImp        = x[end]

                indexImp = findfirst( x -> sol.t[x] == tImp, sol.t)

                Mclose = [ (calcDists(sol.u[i],prob_info.x2,prob_info.rs),i) for i in 2:(length(sol.t)-1) if i != indexImp ]
                sort!(Mclose,by = x -> x[1])
                ts = sol.t[Mclose[1][2]]
                ps = sol.u[Mclose[1][2]]
            end
        end
    end
    return ts, tp, ps, pp
end
function calcOuts(x::Vector{T},sol::ODESolution,par::ProblemParameters,prob_info::ProblemInfo,PROB::Tuple{ODEProblem,ODEProblem}) where {T<:ForwardDiff.Dual}
    ts = sol.t[1] # Instant of close approach to secondary body
    tp = sol.t[end] # Instant of close approach to primary body
    ps = sol.u[1] # State of close approach to secondary body
    pp = sol.u[end] # State of close approach to primary body

    if length(sol.u) > 2
        IMP     = sol.u[2]
    else
        IMP     = sol.u[1]
    end

    if length(x) == 3
        IMP     = sol.u[1]
        if length(sol.t) == 3
            ts = sol.t[2]
            ps = sol.u[2]
        end
    else
        if par.stop
            if length(sol.t) == 3
                ts = sol.t[2]
                ps = sol.u[2]
            end
        else
            if length(sol.t) >= 4
                tImp        = x[end]

                indexImp = findfirst( x -> sol.t[x] == tImp, sol.t)
                IMP  = sol.u[indexImp]

                Mclose = [ (calcDists(sol.u[i],prob_info.x2,prob_info.rs),i) for i in 2:(length(sol.t)-1) if i != indexImp ]
                sort!(Mclose,by = x -> x[1])
                ts = sol.t[Mclose[1][2]]
                ps = sol.u[Mclose[1][2]]
            end
        end
    end
    return ts, tp, ps, pp, IMP
end
"""
Description myObj
"""
function myObj(X::Vector,grad::Vector,par::ProblemParameters)
    if par.stop == 0
        x = X[1:end-1]
    else
        x = X
    end
    #println(x)
    if length(grad) > 0
        if par.stop == 0
            grad[:] = vcat(2x,[0])
        else
            grad[:] = 2x
        end
    end
    return dot(x,x)
end
"""
Description fConstraint
"""
function fConstraint(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem}) where {T<:ForwardDiff.Dual}

    ds, dp, ts, tp  = calc_data(x, par, prob_info, PROB)

    n = par.rest[1]
    if      n == 1
        res = [ dp - par.hgoal ]
    elseif n == 2
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9 ]
    elseif n == 3
        res = [ dp - par.hgoal,
                -tp+par.tspan*0.95 ]
    elseif  n == 4
        res = [ dp - par.hgoal,
                -ds + par.hmmin ]
    elseif n == 8
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9,
                -ds + par.hmmin ]
    end
    return  res
end
function fConstraint(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem,ODEProblem}) where {T<:ForwardDiff.Dual}

    ds, dp, ts, tp, eRet  = calc_data(x, par, prob_info, PROB)

    n = par.rest[1]
    if      n == 1
        res = [ dp - par.hgoal ]
    elseif n == 2
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9 ]
    elseif n == 3
        res = [ dp - par.hgoal,
                -tp+par.tspan*0.95 ]
    elseif  n == 4
        res = [ dp - par.hgoal,
                -ds + par.hmmin ]
    elseif n == 5
        res = [ dp - par.hgoal,
                eRet-1.1*par.hgoal ]
    elseif n == 6
        res = [ dp - par.hgoal,
                -ds+par.hmmin,
                eRet-1.1*par.hgoal ]
    elseif n == 7
        res = [ dp - par.hgoal,
                -eRet,
                eRet-1.1*par.hgoal ]
    elseif n == 8
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9,
                -ds + par.hmmin ]
    end
    return  res
end
function fConstraint(x::Array{Float64,1}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem})

    ds, dp, ts, tp, eRet, ImpF, ImpS, ImpI, tImpF, tImpS, tImpI, fimp  = calc_data(x, par, prob_info, PROB)

    n = par.rest[1]
    if      n == 1
        res = [ dp - par.hgoal ]
    elseif n == 2
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9 ]
    elseif n == 3
        res = [ dp - par.hgoal,
                -tp+par.tspan*0.95 ]
    elseif  n == 4
        res = [ dp - par.hgoal,
                -ds + par.hmmin ]
    elseif n == 5
        res = [ dp - par.hgoal,
                eRet-1.1*par.hgoal ]
    elseif n == 6
        res = [ dp - par.hgoal,
                -ds+par.hmmin,
                eRet-1.1*par.hgoal ]
    elseif n == 7
        res = [ dp - par.hgoal,
                -eRet,
                eRet-1.1*par.hgoal ]
    elseif n == 8
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9,
                -ds + par.hmmin ]
    end
    return  res
end
"""
Description fConstraint
"""
function fAltitude(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem}) where {T<:ForwardDiff.Dual}
    ds, dp, ts, tp  = calc_data(x, par, prob_info, PROB)
    return  dp - par.hgoal
end
function fAltitude(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem,ODEProblem}) where {T<:ForwardDiff.Dual}

    ds, dp, ts, tp, eRet  = calc_data(x, par, prob_info, PROB)

    return  dp - par.hgoal
end
"""
Description fRest
"""
function fRest(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem}) where {T<:ForwardDiff.Dual}

    ds, dp, ts, tp  = calc_data(x, par, prob_info, PROB)

    n = par.rest[1]
    if n == 2
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9 ]
    elseif n == 3
        res = [ dp - par.hgoal,
                -tp+par.tspan*0.95 ]
    elseif  n == 4
        res = [ dp - par.hgoal,
                -ds + par.hmmin ]
    elseif n == 8
        res = [ dp - par.hgoal,
                -dp +  par.hgoal*0.9,
                -ds + par.hmmin ]
    end
    return  res
end
function fRest(x::Vector{T}, par::ProblemParameters,prob_info::ProblemInfo, PROB::Tuple{ODEProblem,ODEProblem}) where {T<:ForwardDiff.Dual}

    ds, dp, ts, tp, eRet  = calc_data(x, par, prob_info, PROB)

    if n == 2
        res = [ -dp +  par.hgoal*0.9 ]
    elseif n == 3
        res = [ -tp+par.tspan*0.95 ]
    elseif  n == 4
        res = [ -ds + par.hmmin ]
    elseif n == 5
        res = [ eRet-1.1*par.hgoal ]
    elseif n == 6
        res = [ -ds+par.hmmin,
                eRet-1.1*par.hgoal ]
    elseif n == 7
        res = [ -eRet,
                eRet-1.1*par.hgoal ]
    elseif n == 8
        res = [ -dp +  par.hgoal*0.9,
                -ds + par.hmmin ]
    end
    return  res
end
"""
Description myConstraint, myConstAlt, myConstRest
"""
function myConstraint(result::Vector,x::Vector,grad::Matrix, par::ProblemParameters,prob_info::ProblemInfo, PROB, results)
    ForwardDiff.jacobian!(results, (x0::Vector) -> fConstraint(x0, par,prob_info, PROB) , x)
    if length(grad) > 0
        grad[:] = DiffResults.jacobian(results)'
    end
    result[:] = DiffResults.value(results)
end
function myConstAlt(result::Vector,x::Vector,grad::Matrix, par::ProblemParameters,prob_info::ProblemInfo, PROB, results)
    ForwardDiff.gradient!(results, (x0::Vector) -> fAltitude(x0, par,prob_info, PROB) , x)
    if length(grad) > 0
        grad[:] = reshape( DiffResults.gradient(results) ,length(x),1)
    end
    result[:] = [DiffResults.value(results)]
end
function myConstRest(result::Vector,x::Vector,grad::Matrix, par::ProblemParameters,prob_info::ProblemInfo, PROB, results)
    ForwardDiff.jacobian!(results, (x0::Vector) -> fRest(x0, par,prob_info, PROB) , x)
    if length(grad) > 0
        grad[:] = DiffResults.jacobian(results)'
    end
    result[:] = DiffResults.value(results)
end
"""
Description fImp
"""
function fImp(p::Vector{<:Real},prob_info::ProblemInfo) # Computes the value of the collinear impulse to perform in p to circularize (adimensional velocity)
    xp, yp, zp, vx, vy, vz = p
    mu  = prob_info.mu  ::Float64
    mup = prob_info.mup ::Float64
    return sqrt((vx-yp)^2+(vy+xp+mu)^2+vz^2) - sqrt( mup / sqrt( (xp+mu)^2 + yp^2 + zp^2 ) )
end
"""
Description fImp
"""
calcDists(y::Vector{<:Real},x::Float64,r::Float64)       = sqrt((y[1]-x)^2+y[2]^2+y[3]^2)-r
"""
Description plotCC
"""
function plotCC(body::Int, n = 16 ::Int, d=prob_info::ProblemInfo) # Function to plot Earth and Moon
    if body == 1 # plot prymary body
        pos             = -d.mu
        rad             = d.rp
        colorSurface    = :blues
    elseif body == 0 # plot secondary
        pos             = 1-d.mu
        rad             = d.rs
        colorSurface    = :Greys
    end
    u = collect(range(0,stop=2,length=n));
    v = collect(range(0,stop=1,length=n));
    x = rad * cospi.(u) * sinpi.(v)' .+ pos;
    y = rad * sinpi.(u) * sinpi.(v)';
    z = rad * repeat(cospi.(v)',outer=[n, 1])
    surface!(x,y,z, color=colorSurface)
end
"""
Description fsObj
"""
function fsObj(x::Vector,grad::Vector, intF)
    if length(grad) > 0
        grad[:] = ForwardDiff.jacobian(intF, x)'
    end
    return intF(x)[1]
end
"""
Description point
"""
function point(i::Int,orbit::Array{Float64,2}) # Uses interpolation to extract the relevant points
    ORB     = interpolate(orbit, BSpline(Quadratic(Reflect(OnCell()))))
    Y(x)    = ORB(x,2)

    if i == 1 # Point of maximum Z coordinate
        imax = argmax(orbit[:,3])
        if (imax==0) || (imax==length(orbit[:,1]))
            ind = imax
        else
            ind = find_zero(Y, imax, Order2())
        end
    elseif i == 4 # Point of minimum Z coordinate
        imax = argmax(-orbit[:,3])
        if (imax==0) || (imax==length(orbit[:,1]))
            ind = imax
        else
            ind = find_zero(Y, imax, Order2())
        end
    elseif i == 6 # Point of most negative Y coordinate
        imax = argmax(-orbit[:,2])
        if (imax==0) || (imax==length(orbit[:,1]))
            ind = imax
        else
            OPT = Opt(:LD_SLSQP, 1)
            #LN_COBYLA
            #LD_SLSQP
            lower_bounds!(OPT, [1.1])
            upper_bounds!(OPT, [200.9])
            ftol_rel!(OPT::Opt, 1e-14::Float64)
            ftol_abs!(OPT::Opt, 1e-14::Float64)
            xtol_rel!(OPT::Opt, 1e-14::Float64)
            xtol_abs!(OPT::Opt, 1e-14::Float64)

            min_objective!(OPT::Opt, (x,grad) -> fsObj(x,grad,Y))
            (optf,optx,ret) = optimize(OPT, [imax])
            ind = optx[1]
        end
    elseif i == 2 # Point of most positive Y coordinate
        imax = argmax(orbit[:,2])
        if (imax==0) || (imax==length(orbit[:,1]))
            ind = imax
        else
            OPT = Opt(:LD_SLSQP, 1)
            lower_bounds!(OPT, [1.1])
            upper_bounds!(OPT, [200.9])
            ftol_rel!(OPT::Opt, 1e-14::Float64)
            ftol_abs!(OPT::Opt, 1e-14::Float64)
            xtol_rel!(OPT::Opt, 1e-14::Float64)
            xtol_abs!(OPT::Opt, 1e-14::Float64)

            max_objective!(OPT::Opt, (x::Vector,grad::Vector) -> fsObj(x,grad,Y))
            (optf,optx,ret) = optimize(OPT, [imax])
            ind = optx[1]

        end
    elseif (i == 3) || (i == 5) # Points of zero Z coordinate
        Z(x)    = ORB(x,3)
        indexes = [i for i in 1:length(orbit[:,1])-1 if orbit[i,3].*orbit[i+1,3] < 0 ]
        if      i == 5
            iz  = [indexes[i] for i in 1:2 if Y(indexes[i]) < 0 ]
            ind = find_zero(Z, iz[1], Order2())
        elseif  i == 3
            iz  = [indexes[i] for i in 1:2 if Y(indexes[i]) > 0 ]
            ind = find_zero(Z, iz[1], Order2())
        end
    end
    return ORB(ind,1:6), ind
end
"""
Description Pset
"""
function Pset(lbounds::Array{Float64,1},ubounds::Array{Float64,1}, x0::Array{Float64,1}, n = 50 ::Int64)::Array{Array{Float64,1},1} # Creates the starting points of the GlobalSearch like MultiStart
    cuts = Int(floor(n^(1/length(x0)))) # number of cuts in the grid
    C = ones(Int, length(x0))*(cuts + 1)
    p = prod(C); i = 1
    while p > n
        C[i]    -= 1
        p       = prod(C)
        i       += 1
        if i > length(x0)
            i = 1
        end
    end

    for i in 1:length(C)
        if C[i] == 1
            C[i] += 1
        end
    end

    rangeVar = [collect(range(lbounds[i],stop = ubounds[i], length = C[i])) for i in 1:length(x0)]
    if length(x0) == 3
        setOfPoints = [ [i, j, k] for i=rangeVar[1] for j=rangeVar[2] for k=rangeVar[3]  ]
    elseif length(x0) == 4
        setOfPoints = [ [i,j,k,l] for i=rangeVar[1] for j=rangeVar[2] for k=rangeVar[3] for l=rangeVar[4]  ]
    elseif length(x0) == 5
        setOfPoints = [ [i,j,k,l,m] for i=rangeVar[1] for j=rangeVar[2] for k=rangeVar[3] for l=rangeVar[4] for m=rangeVar[5]  ]
    elseif length(x0) == 6
        setOfPoints = [ [i,j,k,l,m,n] for i=rangeVar[1] for j=rangeVar[2] for k=rangeVar[3] for l=rangeVar[4] for m=rangeVar[5] for n=rangeVar[6]  ]
    elseif length(x0) == 7
        setOfPoints = [ [i,j,k,l,m,n,p] for i=rangeVar[1] for j=rangeVar[2] for k=rangeVar[3] for l=rangeVar[4] for m=rangeVar[5] for n=rangeVar[6] for p=rangeVar[7]  ]
    else
        println("Wrong input, x0 dimension greater than 7")
        return [[0.0]]
    end

    while length(setOfPoints) < n-1
        point = ubounds.*rand(length(x0))+lbounds.*rand(length(x0))
        setOfPoints = cat(setOfPoints, [point], dims = 1)
    end
    return cat(setOfPoints, [x0], dims = 1)
end
"""
Description setOptProblem
"""
function setOptProblem(par::ProblemParameters,prob_info::ProblemInfo,PROB)

    if par.stop == 1
        lbounds = -ones(length(par.x0))*par.lim
        ubounds =  ones(length(par.x0))*par.lim
    else
        if par.maxt0 < 0
            lbounds = vcat(-ones(length(par.x0)-1)*par.lim, [par.maxt0] )
            ubounds = vcat( ones(length(par.x0)-1)*par.lim, [-1e-3] )
        else
            lbounds = vcat(-ones(length(par.x0)-1)*par.lim, [-1e-3] )
            ubounds = vcat( ones(length(par.x0)-1)*par.lim, [par.maxt0] )
        end
    end

    # Optimization
    Nopt = length(par.x0) # Number of optimization Variables
    alg = :LD_SLSQP
    OPT = Opt(alg, Nopt)
    lower_bounds!(OPT::Opt, lbounds)
    upper_bounds!(OPT::Opt, ubounds)
    min_objective!(OPT::Opt, (X::Vector,grad::Vector) -> myObj(X,grad,par))
    ftol_rel!(OPT::Opt, par.optTol)
    ftol_abs!(OPT::Opt, par.optTol)
    xtol_rel!(OPT::Opt, par.optTol)
    xtol_abs!(OPT::Opt, par.optTol)
    maxeval!(OPT::Opt, par.neval)

    if par.tRest
        # Create DiffResults object
        results = DiffResults.JacobianResult(ones(par.rest[2]),par.x0)
        inequality_constraint!(OPT, (result::Vector,x::Vector,grad::Matrix) -> myConstraint(result,x,grad, par,prob_info, PROB, results), par.cTol*ones(par.rest[2]) )
    else
        if par.rest[2] > 1
            results = DiffResults.JacobianResult(ones(par.rest[2]-1),par.x0)
            inequality_constraint!(OPT, (result::Vector,x::Vector,grad::Matrix) -> myConstRest(result,x,grad, par,prob_info, PROB, results), par.cTol*ones(par.rest[2]-1) )
        end
        Eresults = DiffResults.GradientResult(par.x0)
        equality_constraint!(OPT, (result::Vector,x::Vector,grad::Matrix) -> myConstAlt(result,x,grad, par,prob_info, PROB, Eresults), par.cTol*ones(1) )
    end

    if par.sType
        return OPT, Pset(lbounds,ubounds,par.x0,par.points)
    else
        return OPT, [par.x0]
    end
end
"""
Description createPROB
"""
function createPROB(par::ProblemParameters, prob_info::ProblemInfo)
    call, cRet          = createCallback(par,prob_info)
    # create ODE problem
    if par.retorno
        prob    = ODEProblem((dy,y,p,t) -> dineq!(dy,y,p,t,prob_info.mu), par.y, par.tspan, [0.0, 0.0, 0.0, 0.0, 1.0], callback = call)
        probRet = ODEProblem((dy,y,p,t) -> dineq!(dy,y,p,t,prob_info.mu), par.y, par.tspan, [0.0, 0.0, 0.0, 0.0, 1.0], callback = cRet)
        PROB    = tuple(prob,probRet)
    else
        prob    = ODEProblem((dy,y,p,t) -> dineq!(dy,y,p,t,prob_info.mu), par.y, par.tspan, [0.0, 0.0, 0.0, 0.0, 1.0], callback = call)
        PROB    = tuple(prob)
    end
    return PROB
end
"""
Description createPROB
"""
function search(par::ProblemParameters, prob_info::ProblemInfo)
    # create ODE problem(s)
    PROB = createPROB(par, prob_info)
    #set optimization problem
    OPT, searchPoints   = setOptProblem(par,prob_info,PROB)
    @show OPT
    # Perform search
    sol = map((x) -> myOptimizer(x,OPT),searchPoints)
    return sol
end
"""
Description createCallback
"""
function createCallback(par::ProblemParameters,prob_info::ProblemInfo)
    lX0     = length(par.x0)
    t2      = par.tspan
    detMin  = par.detMin
    stop    = par.stop
    col     = par.col
    mu      = prob_info.mu


    cRet   = ContinuousCallback((y,t,integrator) -> eConditionRet(y,t,integrator,mu),terminate!, nothing, save_positions = (false,true), rootfind=true) # Fee Return stoping

    if lX0 == 3

        if t2 < 0
            cDist   = ContinuousCallback((y,t,integrator) -> eCondition(y,t,integrator,mu),nothing, terminate!, save_positions = (false,true), rootfind=true)
            mDist   = ContinuousCallback((y,t,integrator) -> mCondition(y,t,integrator,mu),nothing, saveEvent!, save_positions = (false,true), rootfind=true)
        else
            cDist   = ContinuousCallback((y,t,integrator) -> eCondition(y,t,integrator,mu),terminate!,nothing, save_positions = (false,true), rootfind=true)
            mDist   = ContinuousCallback((y,t,integrator) -> mCondition(y,t,integrator,mu),saveEvent!, nothing, save_positions = (false,true), rootfind=true)
        end

        if detMin == 1
            call = CallbackSet(cDist,mDist)
        else
            call = CallbackSet(cDist)
        end

    else

        if t2 < 0
            if stop == 1
                cDist   = ContinuousCallback((y,t,integrator) -> gCondition(y,t,integrator,mu),nothing, (integrator) -> impulse!(integrator,col),   save_positions = (false,true), rootfind=true)
            else
                mDist   = ContinuousCallback((y,t,integrator) -> mCondition(y,t,integrator,mu),nothing, saveEvent!, save_positions = (false,true), rootfind=true)
                eDist   = ContinuousCallback((y,t,integrator) -> eCondition(y,t,integrator,mu),nothing, terminate!, save_positions = (false,true), rootfind=true)
            end
        else
            if stop == 1
                cDist   = ContinuousCallback((y,t,integrator) -> gCondition(y,t,integrator,mu), (integrator) -> impulse!(integrator,col),   nothing, save_positions = (false,true), rootfind=true)
            else
                mDist   = ContinuousCallback((y,t,integrator) -> mCondition(y,t,integrator,mu),saveEvent!, nothing, save_positions = (false,true), rootfind=true)
                eDist   = ContinuousCallback((y,t,integrator) -> eCondition(y,t,integrator,mu),terminate!, nothing, save_positions = (false,true), rootfind=true)
            end
        end

        if stop == 1
            call = CallbackSet(cDist)
        else
            iDist   = DiscreteCallback(iCondition, (integrator) -> impulseInt!(integrator,col),save_positions = (false,true) )
            if detMin == 0
                call = CallbackSet(eDist,iDist)
            else
                call = CallbackSet(eDist,iDist,mDist)
            end
        end
    end
    return call, cRet
end
"""
Description conditions and affect!'s
"""
function gCondition(y,t,integrator,mu::Float64)
    xp, yp, zp, vx, vy, vz = y
    #mu      = integrator.p[1]
    # the last parameter indicates where the minimum is to be looked for
    if integrator.p[end] == 1.0  # the minimum is with the secondary body when p[end] is "true"
        return (xp - 1.0 + mu)*vx +  yp*vy + zp*vz
    else
        return (xp + mu)*vx +  yp*vy + zp*vz
    end
end
function iCondition(y,t,integrator)
    integrator.p[4] == t
end
function eCondition(y::Array{Float64,1},t::Float64,integrator,mu::Float64)  #Float version
    xp, yp, zp, vx, vy, vz = y
    #mu      = integrator.p[1]
    tint    = integrator.p[4]
    if xp > 0.6 || t > tint
        return 0.0
    else
        return (xp + mu)*vx +  yp*vy + zp*vz
    end
end
function eCondition(y,t,integrator,mu::Float64)
    xp, yp, zp, vx, vy, vz = y
    #mu      = integrator.p[1]
    tint    = integrator.p[4]
    if xp > 0.6 || t > tint
        return convert(eltype(y),0.0)
    else
        return (xp + mu)*vx +  yp*vy + zp*vz
    end
end
function eConditionRet(y::Array{Float64,1},t::Float64,integrator,mu::Float64) # Float version
    xp, yp, zp, vx, vy, vz = y
    #mu = integrator.p[1]
    if xp > 0.6
        return 0.0
    else
        return (xp + mu)*vx +  yp*vy + zp*vz
    end
end
function eConditionRet(y,t,integrator,mu::Float64)
    xp, yp, zp, vx, vy, vz = y
    #mu = integrator.p[1]
    if xp > 0.6
        return convert(eltype(y),0.0)
    else
        return (xp + mu)*vx +  yp*vy + zp*vz
    end
end
function mCondition(y,t,integrator,mu::Float64)
    xp, yp, zp, vx, vy, vz = y
    #mu = integrator.p[1]
    cond = (xp - 1.0 + mu)*vx +  yp*vy + zp*vz
    return cond
end
function impulse!(integrator,col::Bool)
    imp = integrator.p[1:3]
    if integrator.p[end] == 1.0
        if col
            integrator.u[4:6] = integrator.u[4:6]*(1.0-imp[1]/norm(integrator.u[4:6]))
        else
            integrator.u[4:6] = integrator.u[4:6] - imp
        end
        integrator.p[end] = convert(eltype(integrator.u),0.0) # change the minimum search to the other body
    else
        terminate!(integrator)
    end
end
function saveEvent!(integrator)
    nothing
    #println("Moon encounter at t: ",integrator.t)
end
function impulseInt!(integrator,col::Bool)
    imp = integrator.p[1:3]
    if col
        integrator.u[4:6] = integrator.u[4:6]*(1.0-imp[1]/norm(integrator.u[4:6]))
    else
        integrator.u[4:6] = integrator.u[4:6] - imp
    end
end
"""
Description myOptimizer
"""
function myOptimizer(x::Array{Float64,1}, OPT)
    f, xf, res = optimize(OPT,x)
    println("The optimizer reached the value: ",f)
    println(res)
    return f, xf, res, x
end
"""
Description myOptimizer
"""
function filterGlobalSol(Gsol, par::ProblemParameters,prob_info::ProblemInfo)
    # create callbacks
    PROB = createPROB(par, prob_info)
    #
    #Eresults = DiffResults.GradientResult(par.x0)
    #for i in 1:length(Gsol)
    #    ForwardDiff.gradient!(Eresults, (x0::Vector) -> fAltitude(x0, par,prob_info, PROB) , Gsol[i][2])
    #    @show abs.( [DiffResults.value(Eresults)] ), Gsol[i][1]

        #@show abs.( fConstraint( Gsol[i][2], par,prob_info, PROB)  ), Gsol[i][1]
    #end
    filterSol = [ Gsol[i] for i in 1:length(Gsol) if abs.( fConstraint( Gsol[i][2], par,prob_info, PROB) ) < par.cTol*ones(par.rest[2]) ]

    sfSol = [filterSol[1]] # filtrar por diferencias
    for i in 2:length(filterSol)
        #println(filterSol[i][1])
        dif = 1
        for j in 1:length(sfSol)
            difF = abs(sfSol[j][1] - filterSol[i][1])/maximum([1,sfSol[j][1]])
            difX = norm(sfSol[j][2] - filterSol[i][2])/maximum([1,norm(sfSol[j][2])])
            if (difF < 0.1) && (difX < 0.1)
                dif = 0
            end
        end
        if dif == 1
            #println("new sol")
            sfSol = vcat(sfSol, filterSol[i])
        end
    end
    println("Number of solutions: ",length(sfSol))
    sfSol = sort!(sfSol, by = x -> x[1])
    return filterSol, sfSol
end
"""
Description plotSols
"""
function plotSols(sfSol, par::ProblemParameters, prob_info::ProblemInfo, orbit::Matrix{<:Real},guiP::Bool)
    # create ODE problem
    PROB = createPROB(par, prob_info)

    pyplot()
    ploty = plot(orbit[:,1],orbit[:,2],orbit[:,3], yflip=true, leg=true, lab="Halo Orbit", aspect_ratio=:equal)
    scatter!([-prob_info.mu],[0],[0],c=:blue,label="Earth",markersize=2)
    scatter!([1-prob_info.mu],[0],[0],c=:grey,label="Moon",markersize=1)

    y0 = par.y
    prob = PROB[1]
    for i in 1:length(sfSol)

        x = sfSol[i][2]

        ds, dp, ts, tp, eRet, ImpF, ImpS, ImpI, tImpF, tImpS, tImpI, fimp = calc_data(x, par, prob_info, PROB)

        if par.stop
            if length(x) == 3
                intProb           = remake(prob, u0 = y0 - vcat(zeros(3),x[1:3]) )
            elseif length(x) == 4
                intProb           = remake(prob, u0 = y0 - vcat(zeros(3),x[1:3]) , p = [x[4], 0.0, 0.0, 0.0, par.detMin] )
            else
                intProb           = remake(prob, u0 = y0 - vcat(zeros(3),x[1:3]) , p = [x[4], x[5], x[6], 0.0, par.detMin] )
            end
            tstop = [0.0]
        else
            if length(x) == 5
                intProb           = remake(prob, u0 = y0 - vcat(zeros(3),x[1:3]) , p = [x[4], 0.0, 0.0, x[5], par.detMin] )
            elseif length(x) == 7
                intProb           = remake(prob, u0 = y0 - vcat(zeros(3),x[1:3]) , p = [x[4], x[5], x[6], x[7], par.detMin] )
            end
            tstop = [x[end]]
        end
        if par.odeAlg
            sol         = solve(intProb, DP5(), abstol= par.atol, reltol= par.rtol, tstops =  tstop , dense=true, save_everystep=true)
        else
            sol         = solve(intProb, Vern7(), abstol= par.atol, reltol= par.rtol, tstops = tstop , dense=true, save_everystep=true)
        end
        if length(x) == 3
            deltav = @sprintf("%.2f (m/s)",(fimp + norm(x))*prob_info.r12*prob_info.w)
        elseif par.stop
            deltav = @sprintf("%.2f (m/s)",(fimp + norm(x[1:3]) + norm(x[4:end]))*prob_info.r12*prob_info.w)
        else
            deltav = @sprintf("%.2f (m/s)",(fimp + norm(x[1:3]) + norm(x[4:end-1]))*prob_info.r12*prob_info.w)
        end
        plot!(sol, vars=(1,2,3), leg=true, label=string(deltav), xlim=:auto)
    end
    #xlims!(-0.05,1.1)
    if guiP
        gui()
        return ploty
    end
end
