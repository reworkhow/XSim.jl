"""
    child = get_child(father::Animal,mother::Animal)

* Produce one child from a mating between the **father** and the **mother**.

"""
function get_child(father::Animal, mother::Animal)
    return Animal(father, mother)
end

function random_mate(father::Animal, mother::Animal)
    return mate(father, mother)
end


function embryo_transfer(mother::Animal, father::Animal)
    return mate(mother, father)
end

## List requested functions as pseudo codes below

# function get_progenies(dams::Cohort, sires::Cohort)
# end

# function get_progeny(dam::Animal, sire::Animal)
# end

# function select(cohort::Cohort, n_select::Int64;
#                 criteria::String="phenotypes", is_positive_sel::Bool=true)

#     if n_select > cohort.n
#         @warn "Selection number is capped to the cohort size."
#         n_select = cohort.n
# end

# function select(cohort::Cohort, prop_select::Float16;
#                 criteria::String="phenotypes", is_positive_sel::Bool=true)
#     # get n to select from the input proportion
#     n_select = round(cohort.n * prop_select)

#     return select(cohort, n_select;
#                   criteria=criteria, is_positive_sel=is_positive_sel)
# end
