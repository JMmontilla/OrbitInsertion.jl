using OrbitalInsertion, InterpHalo, BenchmarkTools

# Create the Halo orbit with InterpHalo
ORB = map(x -> intH(x,0.5,2), range(0.0,stop=1.0,length=201) )
orbita = zeros(length(ORB),6)
for i in 1:size(orbita,1)
    orbita[i,:] = ORB[i]
end
# Select one point of the orbit
y0 = point(1,orbita)[1]

# Create the problem relevant magnitudes
sys         = createSystem(EMsystem()) # when EMsystem is used it gives back the properties of the Earth-Moon system
prob_info   = createConstants(sys)  # which are now used to create the problem information object needed for the optimization

firstGuess = [0.05,-0.05,0.1]
par = createPar(prob_info;  tspan = -2.8, hgoal = 300, hmmin = 100 ,
                            sType = "global", points = 30, x0 = firstGuess,
                            maxt0 = -1.0, rest = "basic", lim = 3000,
                            cTol = 1e-5, tRest = "mixed", alg = "lowOrder",
                            atol = 1e-9, rtol = 1e-9, retorno = false, detMin = false, y = y0, optTol = 1e-8, neval = 450)

gsol = search(par,prob_info)
# @time search(par,prob_info) # test the time to compute the global search

fsol, sFsol = filterGlobalSol(gsol, par,prob_info)

plot = plotSols(sFsol, par, prob_info, orbita,true)

# Make local search with the best of the global solutions
X0 = sFsol[2][2]

par = createPar(prob_info;  tspan = -2.8, hgoal = 300, hmmin = 100 ,
                            sType = "local", x0 = X0,
                            maxt0 = -1.0, rest = "basic", lim = 3000,
                            cTol = 1e-10, tRest = "mixed", alg = "lowOrder",
                            atol = 1e-14, rtol = 1e-14, retorno = false, detMin = false, y = y0, optTol = 1e-12, neval = 1500)

finalSolution = search(par,prob_info)
plot = plotSols(finalSolution, par, prob_info, orbita,true)
