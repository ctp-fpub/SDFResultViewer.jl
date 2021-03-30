function fix_text!(ax, default=:small)
    sc = ax.scene[OldAxis]
    if default == :big
        sc[:names, :textsize] = (25, 25, 25)
        sc[:ticks, :textsize][] = (15, 20, 20)
        sc[:frame, :linewidth][] = (180, 180, 180)
    else
        sc[:names, :textsize] = (10, 10, 10)
        sc[:ticks, :textsize][] = (10, 10, 10)
        sc[:frame, :linewidth][] = (180, 180, 180)
    end
end

function default3d_plot!(fig, f::ScalarField, cmap; kwargs...)
    if !haskey(kwargs, :levels)
        kwargs = Iterators.flatten((kwargs, (:levels=>6,)))
    end

    contour(fig[1,1], f,
        alpha = 0.3, transparency = true, colormap=cmap; kwargs...)
end

function default3d_plot!(fig, f::ScalarVariable, cmap; kwargs...)
    size = get(kwargs, :size3d, 35)

    scattervariable(fig[1,1], f, size=size)
end

function section_plot!(fig, f::ScalarField, slice_dir, idx_value; kwargs...)
    dim = dir_to_idx(slice_dir)
    levels = get(kwargs, :levels, 6)
    f_section = @lift(slice(f, dim, $idx_value))
    labels = string.(filter(i->i≠slice_dir, (:x,:y,:z)))

    section_ax, section_plt = contour(fig[1,2], f_section, levels=levels)

    section_ax.xlabel = labels[1]
    section_ax.ylabel = labels[2]
    section_ax.aspect = DataAspect()

    return section_ax, section_plt
end

function section_plot!(ax, f::ScalarVariable, slice_dir, idx_value; kwargs...)
    dim = dir_to_idx(slice_dir)
    ϵ = get(kwargs, :ϵ, 0.2)
    grid = getdomain(f)
    sort!(f, dim)

    f_section = @lift slice(f, dim, $idx_value, ϵ)
    labels = string.(filter(i->i≠slice_dir, (:x,:y,:z)))

    plt = scattervariable!(ax, f_section, size=2)

    ax.xlabel = labels[1]
    ax.ylabel = labels[2]
    ax.aspect = DataAspect()

    return plt
end

function add_legend!(fig, section_plt, ::ScalarField; kwargs...)
    levels = get(kwargs, :levels, 6)
    legend = Colorbar(fig, section_plt,
        label = "Field values",
        width = Relative(3/4),
        height = 10,
        vertical = false,
        tellheight = true,
        colormap = cgrad(:viridis, levels, categorical = true))

    fig[2,1:2] = legend
end

function add_legend!(fig, section_plt, ::ScalarVariable; kwargs...)
    # legend = Colorbar(fig, section_plt,
    #     label = "Variable values",
    #     width = Relative(3/4),
    #     height = 10,
    #     vertical = false,
    #     tellheight = true,
    #     colormap = cgrad(:viridis))

    # fig[2,1:2] = legend
end

function slider(f::ScalarField, slice_dir; kwargs...)
    dim = dir_to_idx(slice_dir)
    grid = getdomain(f)
    slider_domain = grid[dim]
    sl = Slider(axes(slider_domain, 1))
    sl_label = @lift round(slider_domain[$(sl.value)], digits=2)

    return sl, sl_label
end

function slider(f::ScalarVariable, slice_dir; kwargs...)
    dim = dir_to_idx(slice_dir)
    grid = getdomain(f)
    slider_domain = grid[dim]
    sl = Slider(1:2*10^4:length(slider_domain))
    @debug sl.value
    sl_label = @lift round(grid[dim][$(sl.value)], digits=2)

    return sl, sl_label
end

function section_widget(f, slice_dir=:x; kwargs...)
    fig = Figure(resolution=(800, 600))
    f_approx = approximate_field(f)
    grid = getdomain(f_approx)
    dim = dir_to_idx(slice_dir)
    cmap = RGBAf0.(to_colormap(:viridis, 50), 1.0)
    cmap[1:2] .= RGBAf0(0,0,0,0)

    lvl_slider = Slider(1:20)
    # Contour plot
    ax, plt = default3d_plot!(fig, f_approx, cmap; kwargs...)
    fix_text!(ax)

    # Slider
    sl, sl_label = slider(f_approx, slice_dir; kwargs...)

    section_ax = Axis(fig[1,2])
    # Section through the field
    section_plt = section_plot!(section_ax, f_approx, slice_dir, sl.value; kwargs...)

    add_legend!(fig, section_plt, f; kwargs...)

    # Transparent plane

    plane_origin = [minimum(grid)...]
    ws = maximum(grid) .- plane_origin
    ws[dim] = 0

    p = lift(sl.value) do idx
        autolimits!(section_ax)
        update!(fig.scene)

        new_origin = grid[dim][idx]
        plane_origin[dim] = new_origin
        origin = Point3f0(plane_origin...)
        widths = Point3f0(ws...)
        @debug "Plane with origin at $origin"
        FRect3D(origin, widths)
    end

    a = RGBAf0(0,0,0,0)
    c = RGBAf0(0.2, 0.2, 1.0, 1.0)
    img = AbstractPlotting.ImagePattern([c a; a c])
    mesh!(fig[1, 1], p; color = RGBAf0(0.2, 0.2, 1.0, 0.3), transparency=true)

    # Build app
    app = JSServe.App() do
        dom = md"""
        $(fig.scene)
        $(DOM.div("Slice at: ", sl_label, sl))
        """

        DOM.div(JSServe.TailwindCSS, JSServe.Styling, dom)
    end

end
