#!~/Downloads/julia-0.4.5/usr/bin/julia


function readBusData(inputfile)
    local N
    open(inputfile) do filehandle
        key = "BUS DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            N = parse(Int, splitted[1])
        end
    end
    return N
end


function readLoadData(inputfile, N)

    local
        Dn = Array{Set}(N),
        D

    for i = 1:N
        Dn[i] = Set()
    end

    open(inputfile) do filehandle

        key = "LOAD DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            bus_id = parse(Int, splitted[1])
            push!(Dn[bus_id], i) 
            i = i + 1
        end
        D = i - 1

    end
    return Dn, D
end

function readGeneratorData(inputfile, N)

    local
        Gn = Array{Set}(N),                              #Gn[i] = set of generators connected to bus 'i'       
        Gtype = [], 
        G

    open(inputfile) do filehandle
        #READ GENERATOR DATA

        for i = 1:N
            Gn[i] = Set()
        end

        key = "GENERATOR DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            bus_id = parse(Int, splitted[1])
            push!(Gn[bus_id], i) 
            i = i + 1
            push!(Gtype, chomp(splitted[2]))
        end

        G = i - 1
    end
    return Gn, Gtype, G
end

function readGenLookupTable(inputfile)

    local
        Pmin = Dict{AbstractString,Float64}(),           #pu
        Pmax = Dict{AbstractString,Float64}(),           #pu
        RD = Dict{AbstractString,Float64}(),             #pu / min
        RU = Dict{AbstractString,Float64}(),             #pu / min
        DT = Dict{AbstractString,Int}()           #nb of hours
        UT = Dict{AbstractString,Int}()           #nb of hours
        cost = Dict{AbstractString,Float64}()           #€/MWh
        GenOutRate = Dict{AbstractString,Float64}()             #outages/year
        GenOutDuration = Dict{AbstractString,Int}()                  #hours
        GenMaintPeriods = Dict{AbstractString,Int}()            #weeks/year

    open(inputfile) do filehandle
        key = "GENERATOR LOOKUP TABLE"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            gType = chomp(splitted[1])
            push!(Pmin, gType => float(splitted[2]))
            push!(Pmax, gType => float(splitted[3]))
            push!(RD, gType => float(splitted[4]))
            push!(RU, gType => float(splitted[5]))
            push!(DT, gType => float(splitted[6]))
            push!(UT, gType => float(splitted[7]))
            push!(cost, gType => float(splitted[8]))
            push!(GenOutRate, gType => float(splitted[9]))
            push!(GenOutDuration, gType => parse(Int, splitted[10]))
            push!(GenMaintPeriods, gType => parse(Int, splitted[11]))
        end
    end
    return Pmin, Pmax, RD, RU, DT, UT, cost, GenOutRate, GenOutDuration, GenMaintPeriods
end


function readBranchData(inputfile, N)

    local 
        beta = [],
        X = [],
        fmax = [],
        STfmax = [],
        BranchOutRate = [],
        BranchOutDuration = [],
        L
    open(inputfile) do filehandle

        key = "BRANCH DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            temp = zeros(N)
            fr_bus = parse(Int, splitted[1])
            to_bus = parse(Int, splitted[2])
            temp[fr_bus] = 1
            temp[to_bus] = -1
            append!(beta, temp)
            push!(X, float(splitted[3]))
            push!(fmax, float(splitted[4]))
            push!(STfmax, float(splitted[5]))
            push!(BranchOutRate, float(splitted[6]))
            push!(BranchOutDuration, parse(Int, splitted[7]))
        end
        L = length(X)
        beta = reshape(beta, N, L)

    end
    return beta, X, fmax, STfmax, L, BranchOutRate, BranchOutDuration

end

#######Micro-scenario data

function readGenAvailability(inputfile)

    local uGen = []
    open(inputfile) do filehandle

        key = "GENERATOR AVAILABILITY"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(uGen, [parse(Int64,s) for s = splitted])
            i = i + 1
        end
        G = i - 1
        H = round(Int, length(uGen) / G)
        uGen = reshape(uGen, H, G)
        uGen = transpose(uGen)
    end
    return uGen
end

function readBranchAvailability(inputfile)

    local uBranch = []
    open(inputfile) do filehandle

        key = "BRANCH AVAILABILITY"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(uBranch, [parse(Int64,s) for s = splitted])
            i = i + 1
        end
        L = i - 1
        H = round(Int, length(uBranch) / L)
        uBranch = reshape(uBranch, H, L)
        uBranch = transpose(uBranch)
    end
    return uBranch
end


function readMarketOutcome(inputfile)

    local PMCstar = [], uMCstar = []
    open(inputfile) do filehandle

        key = "MARKET GENERATION"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(PMCstar, float(splitted))
            i = i + 1
        end
        G = i - 1
        H = round(Int, length(PMCstar) / G)
        PMCstar = reshape(PMCstar, H, G)
        PMCstar = transpose(PMCstar)

        key = "MARKET ON/OFF"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(uMCstar, [parse(Int64,s) for s = splitted])
            i = i + 1
        end
        G = i - 1
        H = round(Int, length(uMCstar) / G)
        uMCstar = reshape(uMCstar, H, G)
        uMCstar = transpose(uMCstar)
    end
    return PMCstar, uMCstar
end

function readDemandRealization(inputfile)

    local Pload = []     #pu

    open(inputfile) do filehandle

        key = "DEMAND REALIZATION"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(Pload, float(splitted))
            i = i + 1
        end
        D = i - 1
        H = round(Int, length(Pload) / D)
        Pload = reshape(Pload, H, D)
        Pload = transpose(Pload)

    end

    return Pload

end



function readDemandForecast(inputfile)

    local Ploadfc = []     #pu

    open(inputfile) do filehandle

        key = "DEMAND FORECAST"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(Ploadfc, float(splitted))
            i = i + 1
        end
        D = i - 1
        H = round(Int, length(Ploadfc) / D)
        Ploadfc = reshape(Ploadfc, H, D)
        Ploadfc = transpose(Ploadfc)

    end

    return Ploadfc

end


##############
# Data for Micro-scenario generation

function readPeakLoadData(inputfile)

    local
        PCTWeeklyPeakLoad = [],
        PCTDailyPeakLoad = [],
        PCTHourlyPeakLoad = [],
        PCTBusLoad = []

    open(inputfile) do filehandle

        key = "LOAD WEEKLY PEAK DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            push!(PCTWeeklyPeakLoad, float(splitted[1]))
        end

        key = "LOAD DAILY PEAK DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            push!(PCTDailyPeakLoad, float(splitted[1]))
        end

        key = "LOAD HOURLY PEAK DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(PCTHourlyPeakLoad, float(splitted))
        end
        PCTHourlyPeakLoad = reshape(PCTHourlyPeakLoad, 6, 24)
        PCTHourlyPeakLoad = transpose(PCTHourlyPeakLoad)

        key = "BUS LOAD DATA"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            push!(PCTBusLoad, float(splitted[1]))
        end

    end
    return PCTWeeklyPeakLoad, PCTDailyPeakLoad, PCTHourlyPeakLoad, PCTBusLoad
end


#######DA decision data

function readDAdecision(inputfile)

    local PDAstar = [], uDAstar = []
    open(inputfile) do filehandle

        key = "DA GENERATION"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(PDAstar, float(splitted))
            i = i + 1
        end
        G = i - 1
        H = round(Int, length(PDAstar) / G)
        PDAstar = reshape(PDAstar, H, G)
        PDAstar = transpose(PDAstar)

        key = "DA ON/OFF"
        len = length(key)
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= len && line[1:len] == key)
                break
            end
        end
        i = 1
        while !eof(filehandle)
            line = readline(filehandle)
            if(length(line) >= 3 && line[1:3] == "END")
                break
            end
            splitted = split(line,",")
            append!(uDAstar, [parse(Int64,s) for s = splitted])
            i = i + 1
        end
        G = i - 1
        H = round(Int, length(uDAstar) / G)
        uDAstar = reshape(uDAstar, H, G)
        uDAstar = transpose(uDAstar)
    end
    return PDAstar, uDAstar
end
