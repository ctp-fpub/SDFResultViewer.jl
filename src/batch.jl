function batch_plot(sim::EPOCHSimulation, f, folder_name, file_name, args...; overwrite=true, kwargs...)
    if !ispath(joinpath(sim.dir, folder_name))
        mkdir(joinpath(sim.dir, folder_name))
    end

    @time @progress name = "Creating $f plots for $sim.dir" for file in sim
        fn = joinpath(sim.dir, folder_name, file_name(file))
        if isfile(fn)
            @info "File $fn already exists"
            if overwrite
                @debug "Overwriting $fn..."
            else
                @debug "Skipping"
                continue
            end
        end
        plt = f(sim, args...; kwargs...)
        isnothing(plt) && continue
        Plots.savefig(plt, fn)
    end
end

function batch_plot(sims, args...; kwargs...)
    @progress for sim in sims
        batch_plot(sim, args...; kwargs...)
    end
end
