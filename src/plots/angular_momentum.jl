function mean_Lx_plot(sim; species="electron", save=true)
    Lx(file) = compute_L(:x, file, species) |> unit_L
    @time m_Lx = mean(Lx, sim)
    @time m_Lx_plus = mean(Lx, sim, cond=lx->lx>zero(lx))
    @time m_Lx_minus = mean(Lx, sim, cond=lx->lx<zero(lx))
    l = get_parameter(sim, :constant, :m)

    ts = get_time.(sim) .|> unit_t

    plt = Plots.plot(ts, m_Lx,
        title="Average Lx for l = $l",
        xlabel="t",
        ylabel="Lx",
        label="Total Lx")
    Plots.plot!(plt, ts, m_Lx_plus, label="Lx > 0")
    Plots.plot!(plt, ts, m_Lx_minus, label="Lx < 0")

    save && Plots.savefig(plt, "mean_Lx.png")
    plt
end

"""
    Lx_section_plot(file, slice_location; ϵ, species="electron", colormap=:RdYlBu_11, max_size=10^4, label="", save=false)

Plot a section through the Lx of the given `species`, along the direction given by `slice_location`.
The slice has thicknes `ϵ` and the plot is downsampled to `max_size`.
"""
function Lx_section_plot(file, slice_location;
                         ϵ=2e-3unit_l,
                         species="electron",
                         colormap=:RdYlBu_11,
                         max_size=10^4,
                         label="",
                         save=false)

    Lx = uconvert(unit_L, unit_l, compute_L(:x, file, species))
    if all(iszero.(Lx))
        @warn "Lx is 0"
        return nothing
    end
    Lx_section = downsample(slice(Lx, :x, slice_location, ϵ), max_size)
    if isempty(Lx_section)
        @warn "Empty slice"
        return nothing
    end

    cl = max(abs.(extrema(Lx))...)

    bg = Plots.cgrad(colormap, categorical=true)[6]
    t = sprint(show, get_time(file), context=:compact => true)

    plt = Plots.plot(Lx_section;
        xlabel = "y",
        ylabel = "z",
        clims = (-cl,cl),
        cb_title = "Lx",
        title = "Lx at x = $slice_location and t = $t",
        seriescolor = colormap,
        background_inside = bg,
        framestyle = :box,
        label)

    save && Plots.savefig(plt, "Lx_x$(slice_location)_eps$(ϵ)_t$t.png")

    return plt
end
