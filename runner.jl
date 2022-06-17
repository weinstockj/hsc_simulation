include("./simulation.jl")
using CSV

#date = "2021_05_03"
#results = HscSim.main(10000)
#show(results)
#CSV.write("HSC_sim_results_$date.csv", results, delim = "\t")

date = "2021_10_27"
results = HscSim.main(10000, true, 1.20)
show(results)
CSV.write("HSC_sim_results_two_drivers_tet2_strong_$date.csv", results, delim = "\t")

results = HscSim.main(10000, true, 0.80)
show(results)
CSV.write("HSC_sim_results_two_drivers_tet2_weak_$date.csv", results, delim = "\t")
