include("./simulation.jl")
using CSV

date = "2021_05_03"

results = HscSim.main(10000)
show(results)
CSV.write("HSC_sim_results_$date.csv", results, delim = "\t")
