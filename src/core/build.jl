function build(filename ::String;
               Vg       ::Union{Array{Float64 }, Float64}=1.0)
    summary()
end


function Base.summary()
    summary_genome()
    summary_phenome()
end