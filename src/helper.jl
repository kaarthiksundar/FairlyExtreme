# get the output path 
get_output_path() = "./output/"

# get all the case file paths - this is fixed (no cli-arg for the case file and the scenarios)
function get_file_paths() 
    data_path = "./data"
    file_name = "$(data_path)/pglib_opf_case14_ieee.m"
    scenario_file_name = "$(data_path)/scenarios.json"
    valid_scenario_file_name = "$(data_path)/valid_scenarios.csv"
    return (case_file = file_name, 
        scenario_file = scenario_file_name, 
        valid_scenario_file = valid_scenario_file_name)
end 

function get_weight_path(args)
    weighted = args["run_weighted"]
    (weighted == false) && (return "")
    data_path = "./data"
    filename = "$(data_path)/" * args["weight_file"]
    if (!isfile(filename)) 
        @error "file $filename does not exist, quitting..."
        exit()
    end 
    return filename 
end 

# get the indivual loads in the case file 
function get_loads(;weights=Dict()) 
    case, _ = get_file_paths() 
    data = PowerModels.parse_file(case)
    for (i, load) in get(data, "load", [])
        load["weight"] = get(weights, i, 1.0)
    end 
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    return ref[:load] 
end

# get the load served for each individual load based on the solution
function get_load_served(ref, var, m, loads)
    existing_loads = ref[:load] |> keys 
    x_val = JuMP.value.(var[:xd])
    all_loads = loads |> keys 
    load_served = Dict(i => 0.0 for i in all_loads)
    for i in all_loads 
        if !(i in existing_loads)
            continue 
        end 
        load_served[i] = round(x_val[i] * loads[i]["pd"] * 100.0; digits= 4)
    end 
    return load_served
end 

# get the load shed for each individual load based on the solution
function get_load_shed(ref, var, m, loads)
    existing_loads = ref[:load] |> keys 
    x_val = JuMP.value.(var[:xd])
    all_loads = loads |> keys 
    load_shed = Dict(i => 0.0 for i in all_loads)
    isolated_load_shed = 0.0
    for i in all_loads 
        if !(i in existing_loads)
            isolated_load_shed +=  loads[i]["pd"] * 100.0
            continue 
        end 
        load_shed[i] = round((1-x_val[i]) * loads[i]["pd"] * 100.0; digits=4)
    end 
    return load_shed, isolated_load_shed
end

# read all the scenarios in the scenario json file 
function get_scenarios(file)
    scenarios = JSON.parsefile(file)
    return Dict(parse(Int64, k) => val for (k, val) in scenarios)
end 

# read list of valid scenarios from csv 
function get_valid_scenarios(file)
    (!isfile(file)) && (return [])
    return map(x -> Int(x), readdlm(file) |> vec)
end 

# read the weights
function get_weights(file)
    (isempty(file)) && (return Dict{String,Any}())
    return JSON.parsefile(file)
end 

# a function to run the basic formulations to make sure nothing is off
function test_formulations() 
    case, _ = get_file_paths() 
    
    @info "solving DC-OPF for $case"
    run_dc_opf(case)
    
    @info "solving DC-LS for $case"
    run_dc_ls(case; objective_type = :MinSum)
    run_dc_ls(case; objective_type = :MinMax)

    return true
end 

# if min-max causes no load shed, then declare scenario to be a trivial scenario
function check_if_scenario_causes_any_ls(case, scenario, loads, weights)
    _, ref, var, m = run_dc_ls(case, 
        objective_type = :MinMax, 
        p_norm = NaN, 
        scenario = scenario, 
        weights = weights
    )
    ls = get_load_shed(ref, var, m, loads) |> first |> values |> sum
    return !(abs(ls) < 1E-4)
end 

function check_if_scenario_causes_numerical_issues(case, scenario, weights)
    p_norm_values = [2.0, 3.0, 5.0, 10.0]
    for p_norm in p_norm_values
        _, _, _, m = run_dc_ls(case, 
            objective_type = :pNorm, 
            p_norm = p_norm, 
            scenario = scenario, 
            weights = weights
        )
        if termination_status(m) != LOCALLY_SOLVED 
            return false 
        end 
    end     
    return true
end 

# main entry point for executing all the runs based on the command-line arguements
function run(settings; selected_scenarios = [])
    case, scenario_file, valid_scenarios_file = get_file_paths() 
    weight_file = get_weight_path(settings)
    scenarios = get_scenarios(scenario_file)
    valid_scenarios = get_valid_scenarios(valid_scenarios_file)
    weights = get_weights(weight_file)
    loads = get_loads(weights = weights)

    problem_type = settings["problem_type"]
    objective_type = settings["objective_type"]
    p_norm = settings["p_norm"]
    use_constraint = settings["use_constraint"]

    if (use_constraint == false)
        run_model = run_dc_ls
        individual_load_shed = [] 
        filter!(x -> first(x) in valid_scenarios, scenarios)

        for (id, scenario) in scenarios 
            (!isempty(selected_scenarios) && !(id in selected_scenarios)) && (continue)
            valid = check_if_scenario_causes_any_ls(case, scenario, loads, Dict())
            (valid == false) && (continue)
            _, ref, var, m = run_model(case, 
                objective_type = Symbol(objective_type), 
                p_norm = p_norm, 
                scenario = scenario, 
                weights = weights
            )
            load_shed, isolated_load_shed = get_load_shed(ref, var, m, loads)
            sum_load_shed = round(load_shed |> values |> sum; digits=4)  
            total = sum_load_shed + isolated_load_shed
            result = Vector{Any}()
            push!(result, id)
            load_shed_ordered = [load_shed[k] for k in sort(collect(keys(loads)))]
            append!(result, load_shed_ordered)
            append!(result, [sum_load_shed, isolated_load_shed, total])
            push!(individual_load_shed, result)
        end 

        return individual_load_shed
    else 
        run_model = run_dc_ls_with_fairness_constraint 
        individual_load_shed = [] 
        filter!(x -> first(x) in valid_scenarios, scenarios)

        for (id, scenario) in scenarios
            valid = check_if_scenario_causes_any_ls(case, scenario, loads, Dict())
            (valid == false) && (continue)
            _, ref, var, m = run_model(case, 
                scenario = scenario, 
                weights = weights, 
                epsilon = settings["epsilon"]
            )
            (JuMP.termination_status(m) != OPTIMAL) && (continue)
            load_shed, isolated_load_shed = get_load_shed(ref, var, m, loads)
            sum_load_shed = round(load_shed |> values |> sum; digits=4)  
            total = sum_load_shed + isolated_load_shed
            result = Vector{Any}()
            push!(result, id)
            load_shed_ordered = [load_shed[k] for k in sort(collect(keys(loads)))]
            append!(result, load_shed_ordered)
            append!(result, [sum_load_shed, isolated_load_shed, total])
            push!(individual_load_shed, result)
        end 

        return individual_load_shed
    end 
end 
 
# write the results to csv files
function write_results(settings, load_shed)
    problem_type = settings["problem_type"]
    objective_type = settings["objective_type"]
    p_norm = settings["p_norm"]
    use_constraint = settings["use_constraint"]
    epsilon = settings["epsilon"]
    run_weighted = settings["run_weighted"]

    filename = ""
    if use_constraint == false
        filename *= get_output_path() * problem_type * "-"
        if (objective_type in ["MinSum", "MinMax", "MaxLog"])
            filename *= objective_type
        else 
            filename = filename * string(objective_type) * "-" * string(Int(p_norm))
        end 
        (run_weighted) && (filename *= "-w")
        filename *= ".csv"
    else 
        filename *= get_output_path() * "dc-ls-cons-"
        filename *= string(Int(round(epsilon * 100; digits=1)))
        (run_weighted) && (filename *= "-w")
        filename *= ".csv"
    end 

    header = ["scenario_id"]
    append!(header, map(f -> string(f), sort(get_loads() |> keys |> collect)))
    push!(header, "sum_shed")
    push!(header, "isolated_shed")
    push!(header, "total_shed")

    open(filename, "w") do io 
        writedlm(io, [permutedims(header); reduce(hcat, load_shed)'], ',')
    end 
end 

"""
get clean scenarios - removes scenarios that cause numerical issues with IPOPT 
# this can be fixed by outer approximation - not planning to do this for this paper
"""
function get_scenario_ids()
    case, scenario_file = get_file_paths() 
    scenarios = get_scenarios(scenario_file)
    loads = get_loads()
    valid_ids = []
    for (id, scenario) in scenarios 
        valid = check_if_scenario_causes_any_ls(case, scenario, loads, Dict())
        (valid == false) && (continue)
        valid = check_if_scenario_causes_numerical_issues(case, scenario, Dict())
        (valid == false) && (continue)
        push!(valid_ids, id)
    end 
    open("./data/valid_scenarios.csv", "w") do io 
        writedlm(io, valid_ids, ',')
    end
end 

# debug scenario 
function debug_scenario(scenario_id; 
    objective_type = :pNorm, 
    p_norm = 2.0)
    case, scenario_file = get_file_paths() 
    scenarios = get_scenarios(scenario_file)
    loads = get_loads()
    run_model = run_dc_ls

    data, ref, var, m = run_model(case, 
        objective_type = Symbol(objective_type), 
        p_norm = p_norm, 
        scenario = scenarios[scenario_id], 
        debug = true
    )

    load_shed, isolated_load_shed = get_load_shed(ref, var, m, loads)
    total_ls = round(load_shed |> values |> sum; digits=4) + isolated_load_shed
    return (data = data, ref = ref, var = var, m = m, loads = loads, ls = total_ls)
end 

# debug scenario for constrained problem
function debug_scenario(scenario_id, epsilon = 0.0)
    case, scenario_file = get_file_paths() 
    scenarios = get_scenarios(scenario_file)
    loads = get_loads()
    run_model = run_dc_ls_with_fairness_constraint

    data, ref, var, m = run_model(case, 
        epsilon = epsilon,
        scenario = scenarios[scenario_id], 
        debug = true
    )

    load_shed, isolated_load_shed = get_load_shed(ref, var, m, loads)
    total_ls = round(load_shed |> values |> sum; digits=4) + isolated_load_shed
    return (data = data, ref = ref, var = var, m = m, loads = loads, ls = total_ls)
end 
