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

function scalar_field_widget(f::ScalarField, slice_dir=:x; levels=6)
    fig = Figure(resolution=(800, 600))
    f_approx = approximate_field(f)
    cmap = RGBAf0.(to_colormap(:viridis, 50), 1.0)
    cmap[1:2] .= RGBAf0(0,0,0,0)

    lvl_slider = Slider(1:20)
    # Contour plot
    ax, plt = contour(fig[1,1], f_approx,
        levels=levels, alpha = 0.3, transparency = true,
        colormap=cmap)
    fix_text!(ax)

    # Slider
    dim = dir_to_idx(slice_dir)
    domain_range = f_approx.grid[dim]
    sl = Slider(axes(domain_range, 1))
    sl_label = @lift round(f_approx.grid[dim][$(sl.value)], digits=2)
    # sl_format = x->string(round(x, digits=2))
    # sl = labelslider!(fig, "Slice at:", domain_range, format=sl_format)
    # fig[2, 1] = sl.layout

    # Section through the field
    f_section = @lift(slice(f_approx, dim, $(sl.value)))
    labels = string.(filter(i->i≠slice_dir, (:x,:y,:z)))
    section_ax, section_plt = contour(fig[1,2], f_section, levels=levels)

    section_ax.xlabel = labels[1]
    section_ax.ylabel = labels[2]
    section_ax.aspect = DataAspect()

    legend = Colorbar(fig, section_plt,
        label = "Field values",
        width = Relative(3/4),
        height = 10,
        vertical = false,
        tellheight = true,
        colormap = cgrad(:viridis, levels, categorical = true))
    fig[2,1:2] = legend

    # Transparent plane

    plane_origin = [g[1] for g in f_approx.grid]
    ws = [g[end] - g[1] for g in f_approx.grid]
    ws[dim] = 0

    p = lift(sl.value) do idx
        new_origin = f_approx.grid[dim][idx]
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
