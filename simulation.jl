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
end

function HSC()
    return HSC(
        [256],
        [0.],
        5.,
        [false],
        [8.],
        [0]
    )
end

function getBirthRate(cell::HSC, average_fitness::Float64 = 0.)
    fitness = last(cell.fitness) - average_fitness
    if last(cell.count) == 0
        return 0
    else
        return cell.divisionRate * (1 + fitness)
    end
end

function getDeathRate(cell::HSC, average_fitness::Float64 = 0.)
    fitness = last(cell.fitness) - average_fitness
    return cell.divisionRate * (1 - fitness)
end

mutable struct Blood
    n_hsc_cells::Int64
    cells::Vector{HSC}
    n_blood_cells::Vector{Int64}
    driver_mutation_rate::Float64
    passenger_mutation_rate::Float64
    n_drivers::Int64
    average_fitness::Float64
end

function Blood(n_hsc_cells = 100)
    b = Blood(
        n_hsc_cells, # technically these are clones
        Vector{HSC}(undef, n_hsc_cells),
        [n_hsc_cells * 100], # 100 cells per clone
        .0002,
        .005,
        0,
        0.
    )

    for i in range(1, stop = n_hsc_cells)
        b.cells[i] = HSC()
    end

    return b
end

function cycle(b::Blood, dt::Float64)
    n_blood_cells = 0
    average_fitness = 0.
    for c in b.cells
        push!(c.division, last(c.division) + c.divisionRate * dt)
        birthRate = getBirthRate(c, average_fitness) * dt
        deathRate = getDeathRate(c, average_fitness) * dt
        if rand(Bernoulli(birthRate / (birthRate + deathRate)))
            birth(c, birthRate)
        else
            death(c, deathRate)
        end

        if rand(Bernoulli(b.driver_mutation_rate)) && b.n_drivers == 0 # driver mutation rate
            fitness = rand(Beta(4, 16)) # mean of .2
            driverMutate(c, fitness)
            b.n_drivers += 1
        else
            noDriverMutate(c)
        end
        passengerMutate(c, last(c.count) * b.passenger_mutation_rate)
        n_blood_cells += last(c.count)
        average_fitness += last(c.fitness)
    end
    b.average_fitness = average_fitness / b.n_hsc_cells
    push!(b.n_blood_cells, n_blood_cells)
    #printfmt("now have {:d} total blood cells\n", last(b.n_blood_cells))
end

function cycleLifetime(b::Blood, cycles::Int64 = 80, dt::Float64 = 1.)
    for i in 1:cycles
        cycle(b, dt)
    end
end

function noDriverMutate(cell::HSC)
    if !last(cell.hasDriver)
        push!(cell.hasDriver, false)
    else # if already has a driver, can't take it away
        push!(cell.hasDriver, last(cell.hasDriver))
    end
    push!(cell.fitness, last(cell.fitness))
end

function driverMutate(cell::HSC, fitness::Float64 = 0.2)
    push!(cell.hasDriver, true)
    push!(cell.fitness, fitness)
end

function passengerMutate(cell::HSC, rate::Float64)
    if last(cell.count) > 0
        push!(cell.passengers, last(cell.passengers) + rand(Poisson(rate)))
    else
        push!(cell.passengers, last(cell.passengers))
    end
end

function birth(cell::HSC, rate::Float64)
    push!(cell.count, last(cell.count) + rand(Poisson(rate)))
    #printfmt("birth rate {:.2f}; count is now: {:d}\n", rate, last(cell.count))
end

function death(cell::HSC, rate::Float64)
    count = last(cell.count) - rand(Poisson(rate))
    count = count > 0 ? count : 0
    push!(cell.count, count)
    #printfmt("death rate {:.2f}; count is now: {:d}\n", rate, last(cell.count))
end

function summarizeCycles(b::Blood, dt::Float64 = 1.)
    has_driver = false
    driver_index = 0
    driver_birth = 0
    driver_age = 0
    passengers_before_driver = 0
    fitness = 0.
    lastVAF = 0.
    prevVAF = 0.
    deltaVAF = 0.

    # params
    threshold_of_detection = .04
    years_back = 10
    for c in b.cells
        if any(c.fitness .!= 0)
            #println("found driver variant")
            has_driver = true
            driver_index = findall(c.fitness .!= 0)[1] # time of driver event
            driver_birth = driver_index * dt
            driver_age = length(c.division) * dt - driver_birth
            passengers_before_driver = c.passengers[(driver_index - 1)]
            fitness = last(c.fitness)
            lastVAF = last(c.count) / last(b.n_blood_cells)
            if lastVAF < threshold_of_detection
                passengers_before_driver = 0
            end
            prevVAF = c.count[length(c.count) - years_back] / b.n_blood_cells[length(c.count) - years_back]
            deltaVAF = (lastVAF - prevVAF) / float(years_back)
        end
    end

    if has_driver
        return Dict(
            :driver_index => driver_index,
            :driver_age => driver_age,
            :driver_birth => driver_birth,
            :passengers_before_driver => passengers_before_driver,
            :fitness => fitness,
            :lastVAF => lastVAF,
            :prevVAF => prevVAF,
            :deltaVAF => deltaVAF
        )
    else
        return false
    end

end


function main(n_sims::Int64 = 1000)
    results = DataFrame(
        driver_index = Int64[],
        driver_birth = Int64[],
        driver_age = Int64[],
        passengers_before_driver = Int64[],
        fitness = Float64[],
        lastVAF = Float64[],
        prevVAF = Float64[],
        deltaVAF = Float64[]
    )

    for i in 1:n_sims
        if i % 100 == 0
            printfmtln("completed {:d} of {:d} iterations", i, n_sims)
        end
        b = Blood()
        dt = 1.
        cycleLifetime(b, 80, dt)
        r = summarizeCycles(b, dt)
        if typeof(r) != Bool
            push!(results, r)
        end
    end

    return results
end


end
