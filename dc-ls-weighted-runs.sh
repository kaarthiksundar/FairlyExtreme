julia --project=. src/main.jl --problem_type "dc-opf-ls" --objective_type "MinSum" --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.1 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.2 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.3 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.4 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.5 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.6 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.7 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.8 --run_weighted
julia --project=. src/main.jl --use_constraint --epsilon 0.9 --run_weighted