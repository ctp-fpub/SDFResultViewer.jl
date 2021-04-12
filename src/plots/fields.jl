function E_slice_plot(E::ScalarField, dir, slice_location, t;
    issliced=false,
    save=false,
    slice_dir=:x,
    kwargs...)

    if !issliced
        E_slice = slice(E, slice_dir, slice_location)
    else
        E_slice = E
    end

    grid = getdomain(E_slice)
    cl = ustrip(max(abs.(extrema(E_slice))...))
    ts = sprint(show, t, context=:compact => true)
    loc = sprint(show, slice_location, context=:compact => true)
    xlabel, ylabel = filter(x->x≠string(slice_dir), ["x", "y", "z"])

    plt = Plots.heatmap(grid..., E_slice;
        title = "E$dir at t=$ts, x=$loc",
        xlabel, ylabel, colorbar_title = "E$dir",
        seriescolor = :jet1,
        aspect_ratio = 1,
        framestyle = :box,
        clims = (-cl, cl),
        kwargs...)

    if save
        save && Plots.savefig(plt, "E$(dir_slice)_$loc.png")
    end

    return plt
end
