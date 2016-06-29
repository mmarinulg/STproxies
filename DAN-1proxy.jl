#!~/Downloads/julia-0.4.5/usr/bin/julia
tic()
using JuMP
#using GLPKMathProgInterface
using Cbc
using Gurobi
using CPLEX
#using AmplNLWriter

#DEFINE time constants
hpd = 24
dpw = 7
wpy = 52
hpy = hpd * dpw * wpy
dpy = dpw * wpy
Dt = 60                 #minutes
dt = 15                 #minutes


T = []
push!(T, toq())
println("Loading modules took ", T[end])

tic()
#READ ALL NECESSARY DATA
include("./readData.jl")
N = readBusData(string(ARGS[1]))
Dn, D = readLoadData(string(ARGS[1]), N)
Gn, Gtype, G = readGeneratorData(string(ARGS[1]), N)
Pmin, Pmax, RD, RU, UT, DT, cost = readGenLookupTable(string(ARGS[1]))
beta, X, fmax, STfmax, L = readBranchData(string(ARGS[1]), N)
GenSchOut, GenForceOut = readGenAvailability(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), G)
BranchForceOut = readBranchAvailability(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), L)
MarketClearing = readMarketOutcome(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]), G)
Ploadfc = readDemandForecast(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))
push!(T, toq())
println("Reading data took ", T[end])


#CREATE N-1 CONTINGENCY LIST
C = L
a = ones(Float64, L, C) - eye(L, C)

PDAstar = Array(Float64, G, 0)
uDAstar = Array(Int, G, 0)

#DA DECISION
#create empty array indexed by generator, to be filled with tuples (hour, power output)
DAdecision = Array(Any, G)
for g=1:G
    DAdecision[g] = Tuple{Int, Float64}[]
end


#for day=1:dpy
for day=1:1

    tic()

    #first creat vector of generator availability
    uGen = ones(Int, G)
    week = Int(ceil(day/7))
    last_hour = hpd * (day - 1)  #compute last hour in the previous day, to check for forced outages
    for g=1:G
        if week in GenSchOut[g]
            uGen[g] = 0
            continue
        end
        for (start, duration) in GenForceOut[g]
            if start > last_hour break       #if the current outage starts after the considered hour, no need to look at more outages (since GenForceOut is ordered)
            elseif last_hour <= start + duration     #if the generator is out at the considered hour, it is out for the entire day (conservative approach)
                uGen[g] = 0
            end
        end
    end

    #next creat vector of branch availability
    uBranch = ones(Int, L)
    last_hour = hpd * (day - 1)  #compute last hour in the previous day, to check for forced outages
    for l=1:L
        for (start, duration) in BranchForceOut[l]
            if start > last_hour break       #if the current outage starts after the considered hour, no need to look at more outages (since GenForceOut is ordered)
            elseif last_hour <= start + duration     #if the generator is out at the considered hour, it is out for the entire day (conservative approach)
                uBranch[l] = 0
            end
        end
    end
 
    #next creat matrices of market clearing outcome
    uMCstar = zeros(Int, G, hpd)
    PMCstar = zeros(Float64, G, hpd)
    last_hour = hpd * day  #compute last hour of the CURRENT day
    for g=1:G
        for (hour, Pgen) in MarketClearing[g]
            if hour > last_hour break       #no need to look at hours after the current day
            else
                hour_of_day = mod(hour-1,hpd) + 1
                uMCstar[g,hour_of_day] = 1
                PMCstar[g,hour_of_day] = Pgen
            end
        end
    end

    push!(T, toq())
    println("Filling up model data for day ", day, " took ", T[end])

    tic()

    #m = Model(solver = AmplNLSolver("scipampl", ["scip.set"]))
    #m = Model(solver = GurobiSolver())
    #m = Model(solver = GLPKSolverMIP())
    #m = Model(solver = CbcSolver(logLevel=1, threads=16))
    #m = Model()
    m = Model(solver = CplexSolver())


    @variable(m, PDA[1:G,1:hpd])
    @variable(m, uDA[1:G,1:hpd], Bin)
    @variable(m, uMCDA[1:G,1:hpd])
    @variable(m, 0 <= vMCDA[1:G,1:hpd] <= 1)
    @variable(m, 0 <= wMCDA[1:G,1:hpd] <= 1)
    @variable(m, Pup0[1:G,1:hpd] >= 0)
    @variable(m, Pdn0[1:G,1:hpd] >= 0)
    @variable(m, Pupc[1:G,1:C,1:hpd] >= 0)
    @variable(m, Pdnc[1:G,1:C,1:hpd] >= 0)
    @variable(m, f0[1:L,1:hpd])
    @variable(m, th0[1:N,1:hpd])
    @variable(m, STfc[1:L,1:C,1:hpd])
    @variable(m, STthc[1:N,1:C,1:hpd])
    @variable(m, fc[1:L,1:C,1:hpd])
    @variable(m, thc[1:N,1:C,1:hpd])


    @objective(m, Min, sum{sum{cost[Gtype[g]] * (PDA[g,h] + (Pup0[g,h] - Pdn0[g,h])), g=1:G}, h = 1:hpd} )

    #Generator availability and cycling
    @constraint(m, GenAvailability[g=1:G, h = 1:hpd], uDA[g,h] <= uGen[g])
    @constraint(m, DAonoff[g=1:G, h=1:hpd], uDA[g,h] <= 1 - uMCstar[g,h])
    @constraint(m, jointMCDAonoff[g=1:G, h=1:hpd], uMCDA[g,h] == uMCstar[g,h] + uDA[g,h] )
    @constraint(m, jointMCDAupdown1[g=1:G], vMCDA[g,1] - wMCDA[g,1] == uMCDA[g,1] ) #assume all generators are off at the beginning of the day
    @constraint(m, jointMCDAupdown[g=1:G, h=2:hpd], vMCDA[g,h] - wMCDA[g,h] == uMCDA[g,h] - uMCDA[g,h-1] )
    @constraint(m, jointMCDAMinUPTime[g=1:G, h=UT[Gtype[g]]:hpd], sum{vMCDA[g,s], s=h-UT[Gtype[g]]+1:h} <= uMCDA[g,h] )
    @constraint(m, jointMCDAMinDNTime[g=1:G, h=DT[Gtype[g]]:hpd], sum{wMCDA[g,s], s=h-DT[Gtype[g]]+1:h} <= 1 - uMCDA[g,h] ) 

    #Pre-contingency state
    @constraint(m, prePower[n=1:N, h=1:hpd], 
        sum{PMCstar[g,h] + PDA[g,h] + (Pup0[g,h] - Pdn0[g,h]), g in Gn[n]} - 
        sum{beta[n,l] * f0[l,h], l=1:L} == 
        sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
        )
    @constraint(m, preFlow[l=1:L, h=1:hpd], f0[l,h] - uBranch[l] * (1/X[l]) * sum{beta[n,l] * th0[n,h], n=1:N} == 0 )

    @constraint(m, preFlowPosLim[l=1:L, h=1:hpd], f0[l,h] <= fmax[l] )
    @constraint(m, preFlowNegLim[l=1:L, h=1:hpd], -f0[l,h] <= fmax[l] )

    @constraint(m, preGenMinMarket[g=1:G, h=1:hpd], Pdn0[g,h] <= PMCstar[g,h] - uMCstar[g,h] * Pmin[Gtype[g]])
    @constraint(m, preGenMaxMarket[g=1:G, h=1:hpd], Pup0[g,h] <= uMCstar[g,h] * Pmax[Gtype[g]] - PMCstar[g,h])

    @constraint(m, preRampDown[g=1:G, h=1:hpd], Pdn0[g,h] <= 0.5 * Dt * RD[Gtype[g]] )
    @constraint(m, preRampUp[g=1:G, h=1:hpd], Pup0[g,h] <= 0.5 * Dt * RU[Gtype[g]] )

    @constraint(m, DAGenMin[g=1:G, h=1:hpd], PDA[g,h] >= uDA[g,h] * Pmin[Gtype[g]])
    @constraint(m, DAGenMax[g=1:G, h=1:hpd], PDA[g,h] <= uDA[g,h] * Pmax[Gtype[g]])

    @constraint(m, DARampDown[g=1:G, h=2:hpd], PDA[g,h-1] - PDA[g,h] <= Dt * RD[Gtype[g]] )
    @constraint(m, DARampUp[g=1:G, h=2:hpd], PDA[g,h] - PDA[g,h-1] <= Dt * RU[Gtype[g]] )

    #Short-term post-contingency state
    @constraint(m, STpostPower[n=1:N,c=1:C,h=1:hpd], 
        sum{PMCstar[g,h] + PDA[g,h] + (Pup0[g,h] - Pdn0[g,h]), g in Gn[n]} - 
        sum{beta[n,l] * STfc[l,c,h], l=1:L} == 
        sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
        )
    @constraint(m, STpostFlow[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] - uBranch[l] * a[l,c] * (1/X[l]) * sum{beta[n,l] * STthc[n,c,h], n=1:N} == 0 )

    @constraint(m, STpostFlowPosLim[l=1:L,c=1:C,h=1:hpd], STfc[l,c,h] <= STfmax[l] )
    @constraint(m, STpostFlowNegLim[l=1:L,c=1:C,h=1:hpd], -STfc[l,c,h] <= STfmax[l] )

    #Post-contingency state (after the successful application of corrective control)
    @constraint(m, postPower[n=1:N,c=1:C,h=1:hpd], 
        sum{PMCstar[g,h] + PDA[g,h] + (Pup0[g,h] - Pdn0[g,h]) + (Pupc[g,c,h] - Pdnc[g,c,h]), g in Gn[n]} - 
        sum{beta[n,l] * fc[l,c,h], l=1:L} == 
        sum{Ploadfc[d,hpd*(day-1)+h], d in Dn[n]}
        )
    @constraint(m, postFlow[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] - uBranch[l] * a[l,c] * (1/X[l]) * sum{beta[n,l] * thc[n,c,h], n=1:N} == 0 )

    @constraint(m, postFlowPosLim[l=1:L,c=1:C,h=1:hpd], fc[l,c,h] <= fmax[l] )
    @constraint(m, postFlowNegLim[l=1:L,c=1:C,h=1:hpd], -fc[l,c,h] <= fmax[l] )

    @constraint(m, postGenMin[g=1:G,c=1:C,h=1:hpd],
        Pdn0[g,h] + Pdnc[g,c,h] <= PMCstar[g,h] - uMCstar[g,h] * Pmin[Gtype[g]] +
                                    PDA[g,h] - uDA[g,h] * Pmin[Gtype[g]] + Pup0[g,h])
    @constraint(m, postGenMax[g=1:G,c=1:C,h=1:hpd],
        Pup0[g,h] + Pupc[g,c,h] <= uMCstar[g,h] * Pmax[Gtype[g]] - PMCstar[g,h] +
                                    uDA[g,h] * Pmax[Gtype[g]] - PDA[g,h] + Pdn0[g,h])

    @constraint(m, postRampDown[g=1:G,c=1:C,h=1:hpd], Pdnc[g,c,h] <= dt * RD[Gtype[g]] )
    @constraint(m, postRampUp[g=1:G,c=1:C,h=1:hpd], Pupc[g,c,h] <= dt * RU[Gtype[g]] )



    push!(T, toq())
    println("Setting up the model for day ", day, " took ", T[end])

    tic()
    status = solve(m)
    push!(T, toq())
    println("Solving the problem for day ", day, " took ", T[end])

    uDAstar = round(Int, getvalue(uDA))
    PDAstar = getvalue(PDA)

    for g=1:G
        for h=1:hpd
            if uDAstar[g,h] == 1
                hour = hpd * (day - 1) + h
                push!(DAdecision[g], (hour, PDAstar[g,h]))
            end
        end
    end
    
    println("Task ", ARGS[2], ": Objective value = ", getobjectivevalue(m))

end


println("Everything took ", sum(T))


#Write decision to file

outfile = string("./DAdecisions/DAdecision",ARGS[2],".",ARGS[1])
if(isfile(outfile)) rm(outfile) end

open(outfile,"a") do x
    write(x,"DA DECISION\n")
    for g=1:G
        if isempty(DAdecision[g]) continue  end     #skip generators which will remain off
        write(x, string(g), " ")    #write generator index      
        write(x, join(DAdecision[g], " "))
        write(x, "\n")
    end
    write(x,"END\n\n")

end


