function mean_quantity(f, sim)
    tree = maketree(sim.dir => sim.files)

    quantity = FileTrees.load(tree, lazy=true) do file
        println("Loading $(path(file)) on thread $(Threads.threadid())")
        f(file.value)
    end
    exec(mapvalues(mean, quantity))
end

approx_target_size(::ScalarVariable) = 10^6
approx_target_size(::ScalarField) = 160

function approximate_field(f::AbstractPICDataStructure, sz::Int=approx_target_size(f))
    if unit(eltype(f)) ≠ NoUnits
        @debug "Removeing units"
        approximate_field(ustrip(f), sz)
    else
        largest_dim = maximum(size(f))
        factor = largest_dim ÷ sz
        target_size = map(i->i÷factor, size(f))
        @debug "Resizing to size $target_size from $(size(f))"
        subsample(f, target_size...)
    end
end
