using PowerModels 
using JuMP 
using Ipopt 
using JSON
using DelimitedFiles
include("dc-opf.jl")
include("dc-ls.jl")
include("cli-parser.jl")
include("helper.jl")

PowerModels.silence()
 

function cli_parser_run()
    args = get_cli_args() 
    ls = run(args)
    write_results(args, ls)
end 

cli_parser_run()
