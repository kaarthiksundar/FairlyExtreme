using ArgParse

function parse_cli_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--problem_type"
            help = "power flow (dc-ots-ls/dc-opf-ls)"
            arg_type = String 
            default = "dc-opf-ls"
        "--objective_type"
            help = "can be one of MinSum/MinMax/pNorm/MaxLog"
            arg_type = String 
            default = "MinSum"
        "--p_norm", "-p"
            help = "any integer ≥ 2 (do not enter 2, instead enter 2.0) or NaN"
            arg_type = Float64 
            default = NaN
        "--run_weighted" 
            help = "flag to solve weighted case"
            action = :store_true
        "--weight_file"
            help = "file containing weights in data/ folder"
            arg_type = String 
            default = "weights.json"
        "--use_constraint"
            help = "SOCP fairness constraint with min load shed objective; makes --problem_type, --objective_type, --p_norm irrelevant"
            action = :store_true 
        "--epsilon"
            help = "eps value for fairness constraint lies in [0, 1]"
            arg_type = Float64 
            default = 0.0
    end

    return parse_args(s)
end

get_cli_args() = parse_cli_args()