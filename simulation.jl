module HscSim
using Distributions, Random, Formatting, DataFrames, GLM

Random.seed!(1)

mutable struct HSC
    count::Vector{Int64}
    fitness::Vector{Float64}
    divisionRate::Float64 # rate per year
    hasDriver::Vector{Bool}
    division::Vector{Float64} # number of symmetric cell divisions
    passengers::Vector{Int64}
    driverLabel::Vector{Union{Missing, String}}
end

function HSC()
    return HSC(
        [200],
        [0.],
        0.33,
        [false],
        [8.],
        [0],
        [missing]
    )
end

function getBirthRate(cell::HSC, average_fitness::Float64 = 0., dt::Float64 = 1., denom::Float64 = 1.)

    if isnan(average_fitness)
        throw(error("average fitness is nan"))
    end

    fitness = last(cell.fitness) - average_fitness

    if isnan(fitness)
        throw(error("fitness is nan"))
    end

    if last(cell.count) == 0
        return 0
    else
        return dt * last(cell.count) * cell.divisionRate * (1 + fitness) / denom
    end
end

function getDeathRate(cell::HSC, average_fitness::Float64 = 0., dt::Float64 = 1., denom::Float64 = 1.)
    fitness = last(cell.fitness) - average_fitness
    return dt * last(cell.count) * cell.divisionRate * (1 - fitness) / denom
end

mutable struct Blood
    n_hsc_clones::Int64
    clones::Vector{HSC}
    n_blood_cells::Vector{Int64}
    driver_mutation_rate::Float64 # per cell per dt unit of time
    passenger_mutation_rate::Float64 # per cell per dt unit of time
    n_drivers::Int64
    average_fitness::Float64
end

function Blood(n_hsc_clones = 500)
    b = Blood(
        n_hsc_clones, 
        Vector{HSC}(undef, n_hsc_clones),
        [n_hsc_clones * 200],  # 500 clone with 200 cells each to start
        .0002, # not needed
        .006,
        0,
        0.
    )

    for i in range(1, stop = n_hsc_clones)
        b.clones[i] = HSC()
    end

    return b
end

function cycle(b::Blood, dt::Float64, acquire_driver::Bool, driver_alpha::Float64 = 10.0, driver_beta::Float64 = 25.0, driver_label::Union{String, Missing} = missing)

    n_blood_cells = 0
    average_fitness = 0.
    sum_fitness = 0.

    rate_births = 0.

    if acquire_driver # choose among clones that haven't already died out
        clone_idx = Int64[]
        clone_sizes = Int64[]
        idx = 1

        for c in b.clones
            if last(c.count) > 0 && !last(c.hasDriver) # don't give clone multple drivers
                push!(clone_idx, idx)
                push!(clone_sizes, last(c.count))
            end
            idx += 1
        end

        clone_fractions = clone_sizes * 1. / sum(clone_sizes)

        #println("clone_vafs = ")
        #println(clone_vafs)
        #acquiring_clone_idx = clone_idx[rand(1:length(clone_idx))]
        acquiring_clone_idx = clone_idx[rand(Categorical(clone_fractions))]
        #println("Now acquiring driver")
    end

    idx = 1

    for c in b.clones

        push!(c.division, last(c.division) + last(c.count) * c.divisionRate * dt)

        birthRate = getBirthRate(c, average_fitness) * dt
        deathRate = getDeathRate(c, average_fitness) * dt

        if birthRate > 0 && deathRate > 0
            minRate = birthRate / (birthRate + deathRate)

            if minRate < 0. || minRate > 1 || isnan(minRate)
                #throw(DomainError(minRate, "minRate needs to be between 0 and 1"))
                printfmt("birthRate = {:.3f} \n", birthRate)
                printfmt("deathRate = {:.3f} \n", deathRate)
                throw(error("minRate is nan"))
            end

            #printfmt("minRate = {:.3f} \n", minRate)

            if rand(Bernoulli(minRate))
                birth(c, birthRate)
                rate_births += 1.
            else
                death(c, deathRate)
            end
        else
            birth(c, birthRate) # if birthRate and deathRate == 0, then doesn't update clone cell count
        end

        #if rand(Bernoulli(min(1.0, last(c.count) * b.driver_mutation_rate * dt))) && b.n_drivers == 0 # driver mutation rate
        #    fitness = rand(Beta(4, 16)) # mean of .2
        #    driverMutate(c, fitness)
        #    b.n_drivers += 1
        #else
        #    noDriverMutate(c)
        #end

        if acquire_driver && idx == acquiring_clone_idx 
            fitness = rand(Beta(driver_alpha, driver_beta)) # mean of .0025
            #fitness = rand(Beta(5, 1995)) # mean of .0025
            #fitness = rand(Beta(1, 99)) # mean of .01
            driverMutate(c, fitness, driver_label)
            b.n_drivers += 1
        else
            noDriverMutate(c)
        end

        passengerMutate(c, last(c.count) * b.passenger_mutation_rate * dt)
        n_blood_cells += last(c.count)
        sum_fitness += last(c.fitness)
        average_fitness = sum_fitness / idx

        idx += 1
    end
    #b.average_fitness = average_fitness / b.n_hsc_clones
    #printfmt("b.average_fitness = {:.3f} \n", b.average_fitness)
    #printfmt("rate_births = {:.3f} \n", rate_births / length(b.clones))
    push!(b.n_blood_cells, n_blood_cells)
    #printfmt("now have {:d} total blood cells\n", last(b.n_blood_cells))
end

function cycleLifetime(b::Blood, cycles::Int64 = 90, dt::Float64 = 1.)
    cycle_to_acquire_driver = rand(10:cycles)
    for i in 1:cycles
        acquire_driver = cycle_to_acquire_driver == i
        cycle(b, dt, acquire_driver)
    end
end

function cycleLifetimeTwoDrivers(b::Blood, cycles::Int64 = 90, dt::Float64 = 1., relative_fitness_between_drivers = 1.20)
    cycle_to_acquire_dnmt_driver = rand(10:cycles)
    cycle_to_acquire_tet_driver = rand(10:cycles)

    dnmt_driver_alpha = 10.0
    dnmt_driver_beta = 25.0
    tet_relative_fitness = relative_fitness_between_drivers # tet is stronger or weaker than dnmta by this amount
    tet_driver_alpha = dnmt_driver_alpha * tet_relative_fitness
    tet_driver_beta = (dnmt_driver_beta + dnmt_driver_alpha) - tet_driver_alpha
    for i in 1:cycles
        acquire_dnmt_driver = cycle_to_acquire_dnmt_driver == i
        acquire_tet_driver = cycle_to_acquire_tet_driver == i
        if acquire_dnmt_driver
            cycle(b, dt, acquire_dnmt_driver, dnmt_driver_alpha, dnmt_driver_beta, "dnmt")
        elseif acquire_tet_driver
            cycle(b, dt, acquire_tet_driver, tet_driver_alpha, tet_driver_beta, "tet")
        else 
            cycle(b, dt, false)
        end

    end
end

function noDriverMutate(cell::HSC)
    if !last(cell.hasDriver) # no driver
        push!(cell.hasDriver, false)
        push!(cell.driverLabel, missing)
    else # if already has a driver, can't take it away
        push!(cell.hasDriver, last(cell.hasDriver))
        push!(cell.driverLabel, last(cell.driverLabel))
    end
    push!(cell.fitness, last(cell.fitness))
end

function driverMutate(cell::HSC, fitness::Float64 = 0.2, driverLabel::Union{Missing, String} = missing)
    push!(cell.hasDriver, true)
    push!(cell.fitness, fitness)
    push!(cell.driverLabel, driverLabel)
end

function passengerMutate(cell::HSC, rate::Float64)
    if last(cell.count) > 0
        push!(cell.passengers, last(cell.passengers) + rand(Poisson(rate)))
    else
        push!(cell.passengers, last(cell.passengers))
    end
end

function birth(cell::HSC, rate::Float64)
    if rate > 0
        push!(cell.count, last(cell.count) + rand(Poisson(Float64(rate))))
        #printfmt("birth rate {:.2f}; count is now: {:d}\n", rate, last(cell.count))
    else
        push!(cell.count, last(cell.count))
    end
end

function death(cell::HSC, rate::Float64)
    if rate > 0
        count = last(cell.count) - rand(Poisson(Float64(rate)))
    else
        count = last(cell.count)
    end
    count = count > 0 ? count : 0
    push!(cell.count, count)
    #printfmt("death rate {:.2f}; count is now: {:d}\n", rate, last(cell.count))
end

function summarizeCyclesCompetingClones(b::Blood, dt::Float64 = 1., sim_index::Int64 = 1)
    has_driver = false
    driver_index = 0
    driver_birth = 0
    driver_age = 0
    passengers_before_driver = 0
    passengers_before_driver_censored = 0
    fitness = 0.
    lastVAF = 0.
    prevVAF = 0.
    deltaVAF = 0.

    clone_size_to_vaf = 0.5

    clone_count = 0

    coverage = 38
    alt_read_min = 2

    driver_index = Int64[]
    driver_birth = Int64[]
    driver_age = Int64[]
    driver_label = String[]
    clone_count = Int64[]
    fitness = Float64[]
    lastVAF = Float64[]
    passengers_before_driver = Int64[]
    passengers_before_driver_censored = Int64[]

    # params
    years_back = 10
    for c in b.clones
        if any(c.fitness .!= 0)
            #println("found driver variant")
            has_driver = true
            passengers_before_driver_censored_val = 0
            push!(driver_index, findall(c.fitness .!= 0)[1]) # time of driver event
            push!(driver_label, last(c.driverLabel))
            push!(driver_birth, last(driver_index) * dt)
            push!(driver_age, length(c.division) * dt - last(driver_birth))
            push!(passengers_before_driver, c.passengers[(last(driver_index) - 1)])
            push!(fitness, last(c.fitness))
            push!(clone_count, last(c.count))
            push!(lastVAF, clone_size_to_vaf * last(clone_count) / last(b.n_blood_cells))
            if isnan(last(lastVAF))
                printfmt("last(c.count) = {:.3f} \n", last(c.count))
                printfmt("last(b.n_blood_cells) = {:.3f} \n", last(b.n_blood_cells))
            end
            #if lastVAF < threshold_of_detection
            #    passengers_before_driver_censored = 0
            #end

            for p in 1:last(passengers_before_driver)
                passenger_is_observed = rand(Binomial(coverage, last(lastVAF))) >= alt_read_min
                if passenger_is_observed
                    passengers_before_driver_censored_val += 1
                end
            end

            push!(passengers_before_driver_censored, passengers_before_driver_censored_val)

            #prevVAF = clone_size_to_vaf * c.count[length(c.count) - years_back] / b.n_blood_cells[length(c.count) - years_back]
            #deltaVAF = (lastVAF - prevVAF) / float(years_back)
        end
    end

    if has_driver
        return DataFrame(
            simulation_index = repeat([sim_index], length(driver_index)),
            driver_index = driver_index,
            driver_age = driver_age,
            driver_birth = driver_birth,
            driver_label = driver_label,
            passengers_before_driver = passengers_before_driver,
            passengers_before_driver_censored = passengers_before_driver_censored,
            fitness = fitness,
            lastVAF = lastVAF,
            prevVAF = repeat([missing], length(driver_index)), # not implemented yet
            deltaVAF = repeat([missing], length(driver_index)), # not implemented yet
            cloneCount = clone_count,
            nBloodCells = repeat([last(b.n_blood_cells)], length(driver_index))
        )
    else
        return false
    end

end

function summarizeCycles(b::Blood, dt::Float64 = 1.)
    has_driver = false
    driver_index = 0
    driver_birth = 0
    driver_age = 0
    passengers_before_driver = 0
    passengers_before_driver_censored = 0
    fitness = 0.
    lastVAF = 0.
    prevVAF = 0.
    deltaVAF = 0.

    clone_size_to_vaf = 0.5

    clone_count = 0

    coverage = 38
    alt_read_min = 2

    # params
    years_back = 10
    for c in b.clones
        if any(c.fitness .!= 0)
            #println("found driver variant")
            has_driver = true
            passengers_before_driver_censored = 0
            driver_index = findall(c.fitness .!= 0)[1] # time of driver event
            driver_birth = driver_index * dt
            driver_age = length(c.division) * dt - driver_birth
            passengers_before_driver = c.passengers[(driver_index - 1)]
            fitness = last(c.fitness)
            clone_count = last(c.count)
            lastVAF = clone_size_to_vaf * clone_count / last(b.n_blood_cells)
            if isnan(lastVAF)
                printfmt("last(c.count) = {:.3f} \n", last(c.count))
                printfmt("last(b.n_blood_cells) = {:.3f} \n", last(b.n_blood_cells))
            end
            #if lastVAF < threshold_of_detection
            #    passengers_before_driver_censored = 0
            #end

            for p in 1:passengers_before_driver
                passenger_is_observed = rand(Binomial(coverage, lastVAF)) >= alt_read_min
                if passenger_is_observed
                    passengers_before_driver_censored += 1
                end
            end

            prevVAF = clone_size_to_vaf * c.count[length(c.count) - years_back] / b.n_blood_cells[length(c.count) - years_back]
            deltaVAF = (lastVAF - prevVAF) / float(years_back)
        end
    end

    if has_driver
        return Dict(
            :driver_index => driver_index,
            :driver_age => driver_age,
            :driver_birth => driver_birth,
            :driver_label => missing,
            :passengers_before_driver => passengers_before_driver,
            :passengers_before_driver_censored => passengers_before_driver_censored,
            :fitness => fitness,
            :lastVAF => lastVAF,
            :prevVAF => prevVAF,
            :deltaVAF => deltaVAF,
            :cloneCount => clone_count,
            :nBloodCells => last(b.n_blood_cells)
        )
    else
        return false
    end

end


function main(n_sims::Int64 = 1000, two_drivers = false, relative_fitness_between_drivers = 1.20)
    results = DataFrame(
        driver_index = Int64[],
        driver_birth = Int64[],
        driver_age = Int64[],
        driver_label = Union{String, Missing}[],
        passengers_before_driver = Int64[],
        passengers_before_driver_censored = Int64[],
        fitness = Float64[],
        lastVAF = Float64[],
        prevVAF = Union{Float64, Missing}[],
        deltaVAF = Union{Float64, Missing}[],
        cloneCount = Int64[],
        nBloodCells = Int64[]
    )

    df_list = []

    number_of_cycles = 90
    dt = 1.

    for i in 1:n_sims
        if i % 100 == 0
            printfmtln("completed {:d} of {:d} iterations", i, n_sims)
        end
        b = Blood()
        if two_drivers
            cycleLifetimeTwoDrivers(b, number_of_cycles, dt, relative_fitness_between_drivers)
            r = summarizeCyclesCompetingClones(b, dt, i)
            push!(df_list, r)
        else 
            cycleLifetime(b, number_of_cycles, dt)
            r = summarizeCycles(b, dt)
            if typeof(r) != Bool
                push!(results, r)
            end
        end

    end

    #results = reduce(vcat, pmap(x-> b = Blood();cycleLifetime(b, number_of_cycles, dt); summarizeCycles(b, dt), 1:n_sims))

    if two_drivers
        return reduce(vcat, df_list)
    else
        return results
    end
end


end
