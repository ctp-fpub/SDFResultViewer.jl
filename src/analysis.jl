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
        @debug "Removing units"
        approximate_field(ustrip(f), sz)
    else
        largest_dim = maximum(size(f))
        factor = largest_dim ÷ sz
        factor == 0 && return f
        target_size = map(i->i÷factor, size(f))
        @debug "Resizing to size $target_size from $(size(f))"
        downsample(f, target_size...)
    end
end

function initial_cube(sim, l=100u"nm"; species="electron")
    file = sim[1]
    r = file["grid/$species"]

    findall(rᵢ->rᵢ[1] ∈ zero(l)..l &&
            rᵢ[2] ∈ -l/2..l/2 &&
            rᵢ[3] ∈ -l/2..l/2, r)
end
