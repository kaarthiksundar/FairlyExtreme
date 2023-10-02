#### DC Optimal Power Flow ####
using PowerModels
using Ipopt
using JuMP

function run_dc_opf(file_name::AbstractString)
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)

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
    variables = Dict{Symbol,Any}(
        :va => va, 
        :pg => pg,
        :p => p,
        :p_dc => p_dc  
    ) 

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
    end

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])

    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )


    for (i,bus) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) +                  
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==     
            sum(pg[g] for g in ref[:bus_gens][i]) -                 
            sum(load["pd"] for load in bus_loads) -                 
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        p_fr = p[f_idx]                     
        va_fr = va[branch["f_bus"]]         
        va_to = va[branch["t_bus"]]        
        g, b = PowerModels.calc_branch_y(branch)
        @constraint(model, p_fr == -b*(va_fr - va_to))
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])
    end

    for (i,dcline) in ref[:dcline]
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])  
        @constraint(model, 
            (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    optimize!(model)
    println("The solver termination status is $(termination_status(model))")
    cost = objective_value(model)
    println("The cost of generation is $(cost).")

    return (data = data, ref = ref, variables = variables, model = model)
end 
