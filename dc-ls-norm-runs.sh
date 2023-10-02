julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "MinSum"
julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "pNorm" --p_norm 2.0
julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "pNorm" --p_norm 3.0
julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "pNorm" --p_norm 5.0
julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "pNorm" --p_norm 10.0
julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "MinMax"
