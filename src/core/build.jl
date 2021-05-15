function build(filename ::String;
               Vg       ::Union{Array{Float64 }, Float64}=1.0)
    summary()
end


function summary()
    summary_genome()
    summary_phenome()
end