#### DC Equitable Load Shedding Problem ####
using PowerModels
using Ipopt
using JuMP
using CPLEX

function run_dc_ls(file_name::AbstractString;
    objective_type::Symbol=:MinSum, 
    p_norm::Float64=2.0,
    scenario = [], 
    weights = Dict(), 
    debug = false)
    data = PowerModels.parse_file(file_name)
    for (i, load) in get(data, "load", [])
        load["weight"] = get(weights, i, 1.0)
    end 
    for i in scenario 
        data["branch"][string(i)]["br_status"] = 0 
    end 
    propagate_topology_status!(data)
    
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    model = (objective_type == :MinSum || objective_type == :MinMax) ? Model(CPLEX.Optimizer) : Model(Ipopt.Optimizer)
    
    (debug == false) && (
        if (objective_type == :MinSum || objective_type == :MinMax)
            MOI.set(model, MOI.Silent(), true)
        else
            MOI.set(model, MOI.Silent(), true)
        end
    )

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, 
        ref[:gen][i]["pmin"] <= 
        pg[i in keys(ref[:gen])] <= 
        ref[:gen][i]["pmax"]
    )
    @variable(model, 
        -ref[:branch][l]["rate_a"] <= 
        p[(l,i,j) in ref[:arcs_from]] <= 
        ref[:branch][l]["rate_a"]
    )
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
    @variable(model, p_dc[a in ref[:arcs_dc]])
    @variable(model, 0 <= xd[i in keys(ref[:load])] <= 1)
    variables = Dict{Symbol,Any}(
        :va => va, 
        :pg => pg,
        :p => p,
        :p_dc => p_dc, 
        :xd => xd
    ) 

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
    end

    if (objective_type == :MinSum)
        @objective(model, Min, sum((1 - xd[i]) * load["pd"] * load["weight"] for (i, load) in ref[:load]))
    elseif (objective_type == :MinMax)
        @variable(model, z)
        variables[:max_ls] = z 
        @objective(model, Min, z)
        for (i, load) in ref[:load]
            @constraint(model, z >= load["weight"] * load["pd"] * (1 - xd[i]))
        end 
    elseif (objective_type == :MaxLog)
        @NLobjective(model, Max, sum(log(1e-5 + xd[i] * load["pd"] * load["weight"]) for (i, load) in ref[:load]))
    else 
        @NLobjective(model, Min, sum((1 - xd[i])^p_norm * (load["pd"] * load["weight"])^p_norm for (i, load) in ref[:load]))
    end 
   
    for (i, _) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    for (i, _) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) +                  
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==     
            sum(pg[g] for g in ref[:bus_gens][i]) -                 
            sum(xd[load["index"]] * load["pd"] for load in bus_loads) -                 
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
        )
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        p_fr = p[f_idx]                     
        va_fr = va[branch["f_bus"]]         
        va_to = va[branch["t_bus"]]        
        g, b = PowerModels.calc_branch_y(branch)
        @constraint(model, p_fr == -b*(va_fr - va_to))
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])
    end

    for (i, dcline) in ref[:dcline]
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])  
        @constraint(model, 
            (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    optimize!(model)
    # println("The solver termination status is $(termination_status(model))")

    return (data = data, ref = ref, variables = variables, model = model)
end 

function run_dc_ls_with_fairness_constraint(file_name::AbstractString;
    scenario = [], 
    weights = Dict(), 
    debug = false, 
    epsilon::Float64 = 0.0)

    data = PowerModels.parse_file(file_name)
    for (i, load) in get(data, "load", [])
        load["weight"] = get(weights, i, 1.0)
    end 
    for i in scenario 
        data["branch"][string(i)]["br_status"] = 0 
    end 
    propagate_topology_status!(data)
    
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    model = Model(CPLEX.Optimizer)
    
    (debug == false) && (MOI.set(model, MOI.Silent(), true))

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, 
        ref[:gen][i]["pmin"] <= 
        pg[i in keys(ref[:gen])] <= 
        ref[:gen][i]["pmax"]
    )
    @variable(model, 
        -ref[:branch][l]["rate_a"] <= 
        p[(l,i,j) in ref[:arcs_from]] <= 
        ref[:branch][l]["rate_a"]
    )
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
    @variable(model, p_dc[a in ref[:arcs_dc]])
    @variable(model, 0 <= xd[i in keys(ref[:load])] <= 1)
    @variable(model, t >= 0)
    @variable(model, ls[i in keys(ref[:load])] >= 0)
    @variable(model, t_vec[i in keys(ref[:load])] >= 0)
    variables = Dict{Symbol,Any}(
        :va => va, 
        :pg => pg,
        :p => p,
        :p_dc => p_dc, 
        :xd => xd,
        :t => t,
        :ls => :ls
    ) 

    
    x_vec = []
    n = length(ref[:load]) 
    factor = 1 + (sqrt(n) - 1) * epsilon
    for (i, load) in ref[:load]
        w = load["pd"]
        @constraint(model, ls[i] == (1 - xd[i]) * w)
        @constraint(model, t_vec[i] == factor * ls[i])
        push!(x_vec, t_vec[i])
    end 

    @constraint(model, t == sum(ls)) 
    @constraint(model, [t; x_vec] in SecondOrderCone())

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
    end

    @objective(model, Min, sum((1 - xd[i]) * load["pd"] * load["weight"] for (i, load) in ref[:load]))
   
    for (i, _) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    for (i, _) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) +                  
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==     
            sum(pg[g] for g in ref[:bus_gens][i]) -                 
            sum(xd[load["index"]] * load["pd"] for load in bus_loads) -                 
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
        )
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        p_fr = p[f_idx]                     
        va_fr = va[branch["f_bus"]]         
        va_to = va[branch["t_bus"]]        
        g, b = PowerModels.calc_branch_y(branch)
        @constraint(model, p_fr == -b*(va_fr - va_to))
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])
    end

    for (i, dcline) in ref[:dcline]
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])  
        @constraint(model, 
            (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    optimize!(model)
    # println("The solver termination status is $(termination_status(model))")

    return (data = data, ref = ref, variables = variables, model = model)
end 
