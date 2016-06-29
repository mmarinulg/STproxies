#!~/Downloads/julia-0.4.5/usr/bin/julia

using JuMP
using Cbc
using Gurobi
using CPLEX
using Distributions

T = []

global hpd = 24, dpw = 7, wpy = 52,
		hpy = hpd * dpw * wpy,
		hpw = hpd * dpw,
		dpy = dpw * wpy

Dt = 60

include("./readData.jl")

N = readBusData(string(ARGS[1]))
Dn, D = readLoadData(string(ARGS[1]), N)
Gn, Gtype, G = readGeneratorData(string(ARGS[1]), N)
Pmin, Pmax, RD, RU, DT, UT, cost, GenOutRate, GenOutDuration, GenMaintPeriods = readGenLookupTable(string(ARGS[1]))
beta, X, fmax, STfmax, L, BranchOutRate, BranchOutDuration = readBranchData(string(ARGS[1]), N)
PCTWeeklyPeakLoad, PCTDailyPeakLoad, PCTHourlyPeakLoad, PCTBusLoad = readPeakLoadData(string(ARGS[1]))


#First simulate the annual peak

srand(parse(Int, ARGS[2]))	#seed from argument line
U = Uniform(27, 33) 		#uniformily distributed between 27 and 33 pu
AnnualPeakLoad = rand(U)

#########################
#GENERATOR SCHEDULED OUTAGES

#first create outage array of arrays
GenSchOut = Array(Any, G)

#A) GENERATOR SCHEDULED OUTAGES (eHighWay methodology)

#create array of tuples (generator, capacity) and sort it by decreasing capacity
Gcap = [tuple(g, Pmax[Gtype[g]]) for g=1:G]
orderedGcap = sort(collect(Gcap), by=x->x[2], rev=true)

#initialize virtual residual load (using a weekly basis)
VRLoad = (.01) * PCTWeeklyPeakLoad * AnnualPeakLoad

#main loop
for (g,cap) in orderedGcap
	#println(g, ":", nb_maint[Gtype[g]])
	#Compute number of maintenance periods (of 1 week length)
	n = GenMaintPeriods[Gtype[g]]
	#create array that will be filled with the outages (in weeks)
	GenSchOutTemp = Int[]
	#create a temporary copy of the virtual residual load and compute its maximum (cheat)
	tempVRLoad = VRLoad
	max = maximum(tempVRLoad)
	#loop over maintenance periods
	for m = 1:n
		#pick the week with the lowest virtual residual load
		w = indmin(tempVRLoad)
		#schedule the maintenance
		push!(GenSchOutTemp, w)
		#make sure the same week doesn't get picked again
		tempVRLoad[w] = tempVRLoad[w] + max
	end
	GenSchOut[g] = sort(GenSchOutTemp)
	#now increase VRLoad by the capacity of g in the scheduled periods
	VRLoad[GenSchOut[g]] = VRLoad[GenSchOut[g]] + cap

end

#GENERATOR FORCED OUTAGES

#first create outage array of arrays
GenForceOut = Array(Any, G)

for g=1:G
	#Compute hourly outage rate (outages per hour)
	hGenOutRate = GenOutRate[Gtype[g]] / (hpd * 365)
	#Compute hourly outage probability
	hGenOutProb = hGenOutRate * exp(-hGenOutRate)
	#Compute outage expected duration (in hours)
	ExpGenOutDur = GenOutDuration[Gtype[g]]
	#create array that will be filled with the tuples 'outage time, duration' in hours
	GenForceOut[g] = Tuple{Int,Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		week = Int(ceil(day/7))
		#check if generator already in scheduled maintenance for current week
		if week in GenSchOut[g]
			hour = hour + hpw	#go to next week
			continue
		end
		#simulate forced outage
		if rand(Uniform()) <= hGenOutProb
			#simulate outage duration
			GenOutDur = round(Int, rand(Exponential(ExpGenOutDur)))
			#add tuple
			push!(GenForceOut[g], (hour, GenOutDur))
			hour = hour + GenOutDur
			continue
		end
		hour = hour + 1	
	end
end



#BRANCH FORCED OUTAGES

#first create outage matrix in CRS format
BranchForceOut = Array(Any, L)

for l=1:L
	#Compute hourly outage rate (outages per hour)
	hBranchOutRate = BranchOutRate[l] / (hpd * 365)
	#Compute hourly outage probability
	hBranchOutProb = hBranchOutRate * exp(-hBranchOutRate)
	#Compute outage expected duration (in hours)
	ExpBranchOutDur = BranchOutDuration[l]
	#create array that will be filled with the tuples 'outage time, duration' in hours
	BranchForceOut[l] = Tuple{Int, Int}[]
	#start sweeping all hours in the year
	hour = 1
	while hour <= hpy
		day = Int(ceil(hour/24))
		#check if branch already in scheduled maintenance for current day
		
		#if day in BranchSchOut[l]
		#	hour = hour + hpw	#move to next week
		#	continue
		#end
		
		#simulate forced outage
		if rand(Uniform()) <= hBranchOutProb
			#simulate outage duration
			BranchOutDur = round(Int, rand(Normal(ExpBranchOutDur, 1)))
			#add tuple
			push!(BranchForceOut[l], (hour, BranchOutDur))
			hour = hour + BranchOutDur
			continue
		end
		hour = hour + 1	

	end
end





####################
#DEMAND FORECAST

Ploadfc = Array(Float64, D, 0)

for hour = 1:hpy
	#println("Computing demand realization for hour ", hour)
	day = Int(ceil(hour/24))
	week = Int(ceil(day/7))
	hour_of_day = mod(hour-1,hpd) + 1
	day_of_week = mod(day-1,dpw) + 1
	day_type = day_of_week <= 5 ? 1:2		#day_type: 1 = weekday, 2 = weekend
	if week <= 8 || week >= 44
  		week_type = 0
	elseif week >= 18 && week <= 30
  		week_type = 2
	else
  		week_type = 4
	end 									#week_type: 0 = winter, 2 = summer, 4 = spring or fall
	temp = (.01^4) * PCTHourlyPeakLoad[hour_of_day, week_type + day_type] * PCTDailyPeakLoad[day_of_week] * PCTWeeklyPeakLoad[week] * AnnualPeakLoad * PCTBusLoad
	Ploadfc = hcat(Ploadfc, temp)
end

#################
#DEMAND REALISATION
Pload = Array(Float64, D, 0)

#Laurine's method
sigma_glob = .03
sigma_loc = 0.1
BusLoad = .01 * PCTBusLoad
sum_squares = sum(BusLoad .* BusLoad)

alpha_d = Float16(sqrt((sigma_glob^2 - sum_squares * sigma_loc^2) / (1 - sum_squares)))
beta_d = Float16(sqrt((sigma_glob^2 - sigma_loc^2) / (sum_squares - 1)))

N_alpha = Normal(0, alpha_d)
N_beta = Normal(0, beta_d)

for hour = 1:hpy
	eps_alpha = rand(N_alpha)
	eps_beta = rand(N_beta, D)
	temp = (1 + eps_alpha + eps_beta) .* Ploadfc[:,hour]
	Pload = hcat(Pload, temp)
end


#########################
#MARKET CLEARING OUTCOME

#first create empty array indexed by generator, to be filled with tuples (hour, power output)
MarketClearing = Array(Any, G)
for g=1:G
	MarketClearing[g] = Tuple{Int, Float64}[]
end


#for day=1:dpy
for day=1:1

    tic()

    #first creat vector of generator availability
    uGen = ones(Int, G)
	week = Int(ceil(day/7))
	hour = hpd * (day - 1)	#compute last hour in the previous day, to check for forced outages
    for g=1:G
    	if week in GenSchOut[g]
    		uGen[g] = 0
    		continue
    	end
    	for (start, duration) in GenForceOut[g]
    		if start > hour break		#if the current outage starts after the considered hour, no need to look at more outages (since GenForceOut is ordered)
    		elseif hour <= start + duration 	#if the generator is out at the considered hour, it is out for the entire day (conservative approach)
    			uGen[g] = 0
    		end
    	end
    end
	
	#m = Model(solver = GurobiSolver(OutputFlag=0))
    #m = Model(solver = CbcSolver(logLevel=0, threads=16))
	m = Model(solver = CplexSolver(CPX_PARAM_MIPDISPLAY=0))
	#m = Model(solver = CplexSolver())

	@variable(m, PMC[1:G,1:hpd])
	@variable(m, uMC[1:G,1:hpd], Bin)
	@variable(m, 0 <= vMC[1:G,1:hpd] <= 1)
	@variable(m, 0 <= wMC[1:G,1:hpd] <= 1)

	@objective(m, Min, sum{sum{cost[Gtype[g]] * PMC[g,h], g=1:G}, h = 1:hpd} )

    #Generator availability and cycling
	@constraint(m, GenAvailability[g=1:G, h = 1:hpd], uMC[g,h] <= uGen[g])
	@constraint(m, updown1[g=1:G], vMC[g,1] - wMC[g,1] == uMC[g,1] ) #assume all generators are off at the beginning of the day
	@constraint(m, updown[g=1:G, h = 2:hpd], vMC[g,h] - wMC[g,h] == uMC[g,h] - uMC[g,h-1] )
	@constraint(m, MinUPTime[g=1:G, h = UT[Gtype[g]]:hpd], sum{vMC[g,s], s=h-UT[Gtype[g]]+1:h} <= uMC[g,h] )
	@constraint(m, MinDNTime[g=1:G, h = DT[Gtype[g]]:hpd], sum{wMC[g,s], s=h-DT[Gtype[g]]+1:h} <= 1 - uMC[g,h] )

	#Network balance and limits
	@constraint(m, NetworkBalance[h=1:hpd], sum{PMC[g,h], g=1:G} ==
		sum{Ploadfc[d, hpd * (day - 1) + h], d=1:D} )

	@constraint(m, MinPLim[g=1:G, h = 1:hpd], PMC[g,h] >= uMC[g,h] * Pmin[Gtype[g]] )
	@constraint(m, MaxPLim[g=1:G, h = 1:hpd], PMC[g,h] <= uMC[g,h] * Pmax[Gtype[g]] )

	@constraint(m, RampDown[g=1:G, h=2:hpd], PMC[g,h-1] - PMC[g,h] <= Dt * RD[Gtype[g]] )
	@constraint(m, RampUp[g=1:G, h=2:hpd], PMC[g,h] - PMC[g,h-1] <= Dt * RU[Gtype[g]] )

	status = solve(m)

	uMCstar = round(Int, getvalue(uMC))
	PMCstar = getvalue(PMC)

	for g=1:G
		for h=1:hpd
			if uMCstar[g,h] == 1
				hour = hpd * (day - 1) + h
				push!(MarketClearing[g], (hour, PMCstar[g,h]))
			end
		end
	end

    push!(T, toq())

    println("Task ", ARGS[2], ": Computing market outcome for day ", day, " took ", T[end])

end


#Write data to file

outfile = string("./micro-scenarios/micro",ARGS[2],".",ARGS[1])
if(isfile(outfile)) rm(outfile) end

open(outfile,"a") do x
	write(x,"GENERATOR SCHEDULED OUTAGES\n")
	for g=1:G
		if isempty(GenSchOut[g]) continue end	#skip generators with no outages
		write(x, string(g), " ")	#write generator index	
		write(x, join(GenSchOut[g], " "))
		write(x, "\n")
	end
	write(x,"END\n\n")
	write(x,"GENERATOR FORCED OUTAGES\n")
	for g=1:G
		if isempty(GenForceOut[g]) continue	end 	#skip generators with no outages
		write(x, string(g), " ")	#write generator index	
		write(x, join(GenForceOut[g], " "))
		write(x, "\n")
	end	
	write(x,"END\n\n")
	write(x,"BRANCH FORCED OUTAGES\n")
	for l=1:L
		if isempty(BranchForceOut[l]) continue end 	#skip branches with no outages
		write(x, string(l), " ")	#write branch index	
		write(x, join(BranchForceOut[l], " "))
		write(x, "\n")
	end	
	write(x,"END\n\n")
	write(x,"DEMAND FORECAST\n")
	writedlm(x, Ploadfc, " ")
	write(x,"END\n\n")
	write(x,"DEMAND REALIZATION\n")
	writedlm(x, Pload, " ")
	write(x,"END\n\n")
	
	write(x,"END\n\n")
	write(x,"MARKET CLEARING OUTCOME\n")
	for g=1:G
		if isempty(MarketClearing[g]) continue	end 	#skip generators which will remain off
		write(x, string(g), " ")	#write generator index		
		write(x, join(MarketClearing[g], " "))
		write(x, "\n")
	end
	write(x,"END\n\n")
	
end


