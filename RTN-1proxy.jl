#!~/Downloads/julia-0.4.5/usr/bin/julia
tic()
using JuMP
#using GLPKMathProgInterface
#using Clp
using Gurobi

#constants
global hpd = 24, dpw = 7, wpy = 52,
        hpy = hpd * dpw * wpy,
        dpy = dpw * wpy

#More constants
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
Pmin, Pmax, RD, RU, DT, UT, cost = readGenLookupTable(string(ARGS[1]))
beta, X, fmax, STfmax, L = readBranchData(string(ARGS[1]), N)
PMCstar, uMCstar = readMarketOutcome(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))
PDAstar, uDAstar = readDAdecision(string("./DAdecisions/DAdecision",ARGS[2],".",ARGS[1]))
Pload = readDemandRealization(string("./micro-scenarios/micro",ARGS[2],".",ARGS[1]))

#COMBINE MARKET AND DA
PMCDAstar = PMCstar + PDAstar       #PMCDAstar = PMCstar + PDAstar (uncomment when I have all DA results)
uMCDAstar = uMCstar + uDAstar       #uMCDAstar = uMCstar+ uDAstar (uncomment when I have all DA results)


for g = 1:G
    print(g, "\t")
    print(DT[Gtype[g]], "\t")
    print(UT[Gtype[g]], "\t")
    for h = 1:hpd
        print(uMCstar[g,h], " ")
    end
    println()
    print(g, "\t")
    print(DT[Gtype[g]], "\t")
    print(UT[Gtype[g]], "\t")
    for h = 1:hpd
        print(uMCDAstar[g,h], " ")
    end
    println()
    println("--------------------------------------------------------------")
end



#CREATE N-1 CONTINGENCY LIST
C = L
a = ones(Float64, L, C) - eye(L, C)

#DEFINE time constants
Dt = 60                 #minutes
dt = 15                 #minutes

push!(T, toq())
println("Reading and preparing data took ", T[end])

for h = 1:2*hpd

    tic()
    m = Model(solver = GurobiSolver(OutputFlag=0))
    #m = Model(solver = GLPKSolverLP())
    #m = Model(solver = ClpSolver())
    #m = Model()

    @variable(m, Pup0[1:G] >= 0)
    @variable(m, Pdn0[1:G] >= 0)
    @variable(m, Pupc[1:G,1:C] >= 0)
    @variable(m, Pdnc[1:G,1:C] >= 0)
    @variable(m, f0[1:L])
    @variable(m, th0[1:N])
    @variable(m, STfc[1:L,1:C])
    @variable(m, STthc[1:N,1:C])
    @variable(m, fc[1:L,1:C])
    @variable(m, thc[1:N,1:C])

    @objective(m, Min, sum{cost[Gtype[g]] * (Pup0[g] - Pdn0[g]), g=1:G} )

    #Pre-contingency state
    @constraint(m, prePower[n=1:N], 
        sum{PMCDAstar[g,h] + (Pup0[g] - Pdn0[g]), g in Gn[n]} - 
        sum{beta[n,l] * f0[l], l=1:L} == 
        sum{Pload[d,h], d in Dn[n]}
        )
    @constraint(m, preFlow[l=1:L], f0[l] - (1/X[l]) * sum{beta[n,l] * th0[n], n=1:N} == 0 )

    @constraint(m, preFlowPosLim[l=1:L], f0[l] <= fmax[l] )
    @constraint(m, preFlowNegLim[l=1:L], -f0[l] <= fmax[l] )

    @constraint(m, preGenMin[g=1:G], Pdn0[g] <= PMCDAstar[g,h] - uMCDAstar[g,h] * Pmin[Gtype[g]])
    @constraint(m, preGenMax[g=1:G], Pup0[g] <= uMCDAstar[g,h] * Pmax[Gtype[g]] - PMCDAstar[g,h])

    @constraint(m, preRampDown[g=1:G], Pdn0[g] <= 0.5 * Dt * RD[Gtype[g]] )
    @constraint(m, preRampUp[g=1:G], Pup0[g] <= 0.5 * Dt * RU[Gtype[g]] )

    #Short-term post-contingency state
    @constraint(m, STpostPower[n=1:N,c=1:C], 
        sum{PMCDAstar[g,h] + (Pup0[g] - Pdn0[g]), g in Gn[n]} - 
        sum{beta[n,l] * STfc[l,c], l=1:L} == 
        sum{Pload[d,h], d in Dn[n]}
        )
    @constraint(m, STpostFlow[l=1:L,c=1:C], STfc[l,c] - a[l,c] * (1/X[l]) * sum{beta[n,l] * STthc[n,c], n=1:N} == 0 )

    @constraint(m, STpostFlowPosLim[l=1:L,c=1:C], STfc[l,c] <= STfmax[l] )
    @constraint(m, STpostFlowNegLim[l=1:L,c=1:C], -STfc[l,c] <= STfmax[l] )

    #Post-contingency state (after the successful application of corrective control)
    @constraint(m, postPower[n=1:N,c=1:C], 
        sum{PMCDAstar[g,h] + (Pup0[g] - Pdn0[g]) + (Pupc[g,c] - Pdnc[g,c]), g in Gn[n]} - 
        sum{beta[n,l] * fc[l,c], l=1:L} == 
        sum{Pload[d,h], d in Dn[n]}
        )
    @constraint(m, postFlow[l=1:L,c=1:C], fc[l,c] - a[l,c] * (1/X[l]) * sum{beta[n,l] * thc[n,c], n=1:N} == 0 )

    @constraint(m, postFlowPosLim[l=1:L,c=1:C], fc[l,c] <= fmax[l] )
    @constraint(m, postFlowNegLim[l=1:L,c=1:C], -fc[l,c] <= fmax[l] )

    @constraint(m, postGenMin[g=1:G,c=1:C], Pdn0[g] + Pdnc[g,c] <= PMCDAstar[g,h] - uMCDAstar[g,h] * Pmin[Gtype[g]] + Pup0[g])
    @constraint(m, postGenMax[g=1:G,c=1:C], Pup0[g] + Pupc[g,c] <= uMCDAstar[g,h] * Pmax[Gtype[g]] - PMCDAstar[g,h] + Pdn0[g])

    @constraint(m, postRampDown[g=1:G,c=1:C], Pdnc[g,c] <= dt * RD[Gtype[g]] )
    @constraint(m, postRampUp[g=1:G,c=1:C], Pupc[g,c] <= dt * RU[Gtype[g]] )



    push!(T, toq())
    println("Setting up the model for hour ", h, " took " , T[end])

    #print(m)
    tic()
    status = solve(m)
    push!(T, toq())
    println("Solving the problem for hour ", h, " took " , T[end])

    println("Task ", ARGS[2], ": Objective value = ", getobjectivevalue(m))

end

println("Everything took ", sum(T))


