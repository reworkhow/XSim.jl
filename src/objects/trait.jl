mutable struct Trait
    phenotypic  ::Float64
    genetic     ::Float64
    estimated   ::Float64

    function Trait()
        return new(0.0, 0.0, 0.0)
    end
end