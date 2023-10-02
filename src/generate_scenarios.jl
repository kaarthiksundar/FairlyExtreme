using JSON 
using PowerModels 
using Combinatorics

PowerModels.silence()
data_path = "./data"
file_name = "$(data_path)/pglib_opf_case14_ieee.m"

data = PowerModels.parse_file(file_name)
ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

scenarios = Dict{String,Any}() 

n5 = combinations(ref[:branch] |> keys |> collect, 5) |> collect 

for scenario in n5 
    i = length(scenarios) + 1
    scenarios[string(i)] = scenario 
end 

open("./data/scenarios.json","w") do f
    write(f, JSON.json(scenarios, 2)) 
end