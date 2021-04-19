function approx_plot(f;
    fig = Figure(),
    figurepos = fig[1,1],
    xlabel = "",
    ylabel = "",
    title = "",
    downsample_size = nothing,
    cbar_label="",
    kwargs...)

    ax = Axis(figurepos[1,1];
        xlabel, ylabel,
        title,
        aspect = DataAspect()
    )

    if isnothing(downsample_size)
        f_approx = downsample_approx(f)
    else
        f_approx = downsample_approx(f, downsample_size)
    end
    @debug "Downsampled to $(size(f_approx))"

    plt = fieldplot!(ax, f_approx; kwargs...)

    Colorbar(figurepos[1,2], plt; width=20, label=cbar_label)

    # TODO: Implement save
    # if save
    #     save && Plots.savefig(plt, "E$(dir_slice)_$loc.png")
    # end

    return fig
end

get_unit_label(f) = " (" * unitname(f) * ")"

function compare_slice_with_analytic(numeric_f, f, laser, profile, file;
    fig = Figure(resolution = (800, 700)),
    component=1,
    slice_location=0.0unit_l,
    slice_dir=:x,
    fieldname="",
    f_unit)

    fx = get_parameter(file, :constant, :f_x)*u"m" |> unit_l
    t = get_time(file) |> unit_t
    s = analytic_laser(laser, profile, file)

    grid = getdomain(numeric_f)

    analytic_f = apply_analytic(f, grid, t, s;
        component,
        propagation_dir=:x,
        x₀=fx) |> f_unit

    analytic_slice, loc = dynamic_slice(analytic_f, slice_dir, slice_location)
    numeric_slice, loc = dynamic_slice(numeric_f, slice_dir, slice_location)

    units = get_unit_label(grid)
    xlabel, ylabel = filter(x->x≠string(slice_dir),
        ["x", "y", "z"] .* units)
    ts = sprint(show, t, context=:compact => true)
    loc_str = @lift sprint(show,
        $loc, context=:compact => true) * " " * unitname(slice_location)
    title = @lift fieldname * " slice at t = " * ts * " and " *
        string(slice_dir) * " = " * $loc_str

    ax1 = Axis(fig[1,1][1,1];
        xlabel, ylabel,
        title = "Analytic",
        aspect = DataAspect()
    )
    ax2 = Axis(fig[1,2][1,1];
        xlabel, ylabel,
        title = "Numeric",
        aspect = DataAspect()
    )

    plt1 = fieldplot!(ax1, analytic_slice)
    plt2 = fieldplot!(ax2, numeric_slice)

    Colorbar(fig[1,1][2,1], plt1;
        width = Relative(3/4),
        height = 20,
        vertical = false, flipaxis = false,
        tellheight = true,
        label = fieldname * get_unit_label(numeric_f),
    )
    Colorbar(fig[1,2][2,1], plt2;
        width = Relative(3/4),
        height = 20,
        vertical = false, flipaxis = false,
        tellheight = true,
        label = fieldname * get_unit_label(numeric_f),
    )

    Label(fig[0, :], title, textsize = 30, color = :black)

    return fig
end

function compare_E_slice_with_analytic(laser, profile, file, dir;
    slice_location=0.0unit_l,
    slice_dir=:x)

    E = uconvert(unit_E, unit_l, file["e$dir"])

    compare_slice_with_analytic(E, LaserTypes.E, laser, profile, file;
        component=dir_to_idx(dir),
        slice_location,
        slice_dir,
        f_unit=unit_E,
        fieldname="E$dir"
    )
end

function compare_S_slice_with_analytic(laser, profile, file, dir;
    slice_location=0.0unit_l,
    slice_dir=:x)

    Ex, Ey, Ez = uconvert.(unit_E, unit_l, file[:ex, :ey, :ez])
    Bx, By, Bz = uconvert.(unit_B, unit_l, file[:bx, :by, :bz])

    E = build_vector((Ex, Ey, Ez), (:x, :y, :z))
    B = build_vector((Bx, By, Bz), (:x, :y, :z))
    S = ustrip(E × B) / ustrip(μ_0)

    compare_slice_with_analytic(S, LaserTypes.S, laser, profile, file;
        component=dir_to_idx(dir),
        slice_location=ustrip(slice_location),
        slice_dir,
        f_unit=unit_E,
    )
end

function poynting_summary(file, laser, profile, slice_location)
    t = get_time(file) |> unit_t

    # get data
    fx = get_parameter(file, :constant, :f_x)*u"m" |> unit_l
    Ex, Ey, Ez, E, Bx, By, Bz, B = read_EM_fields(file)

    Ex_grid = getdomain(Ex)
    Bx_grid = getdomain(Bx)

    s = analytic_laser(laser, profile, file)

    analytic_E = get_analytic_E(Ex_grid, t, s, x₀=fx)
    analytic_B = get_analytic_B(Bx_grid, t, s, x₀=fx)
    analytic_S = uconvert(unit_S, analytic_E × analytic_B / μ_0)

    S = uconvert(unit_S, E × B / μ_0)

    S_approx = downsample_approx(S, 15)
    analytic_S_approx = downsample_approx(analytic_S, 15)
    S_slice = downsample_approx(
        slice(S, :x, slice_location),
        50
    )
    analytic_S_slice = downsample_approx(
        slice(analytic_S, :x, slice_location),
        50
    )

    fig = Figure(resolution=(800,800))
    # Layout
    lscene1 = LScene(fig[1, 1], scenekw = (camera = cam3d!, raw = false))
    lscene2 = LScene(fig[1, 2], scenekw = (camera = cam3d!, raw = false))
    ax1 = Axis(fig[2,1][1,1],
        aspect = DataAspect(),
        title="Numeric",
        xlabel="y ($unit_l)",
        ylabel="z ($unit_l)"
    )
    ax2 = Axis(fig[2,2][1,1],
        aspect = DataAspect(),
        title="Analytic",
        xlabel="y ($unit_l)",
        ylabel="z ($unit_l)"
    )
    # Plots
    # 3D
    arrowsize_factor = 0.02
    lengthscale_factor = 10
    fieldplot!(lscene1, S_approx;
        arrowsize_factor,
        lengthscale_factor
    )
    fieldplot!(lscene2, analytic_S_approx;
        arrowsize_factor,
        lengthscale_factor
    )
    # 2D
    arrowsize_factor = 0.15
    lengthscale_factor = 6
    plt1 = fieldplot!(ax1, S_slice;
        arrowsize_factor,
        lengthscale_factor
    )
    plt2 = fieldplot!(ax2, analytic_S_slice;
        arrowsize_factor,
        lengthscale_factor
    )
    # Colorbars
    Colorbar(fig[2,1][2,1], plt1,
        width = Relative(3/4), height = 15,
        vertical = false, flipaxis = false,
        label = "S ($unit_S)",
    )
    Colorbar(fig[2,2][2,1], plt2,
        width = Relative(3/4), height = 15,
        vertical = false, flipaxis = false,
        label = "S ($unit_S)",
    )

    WGLMakie.trim!(fig.layout)
    fig
end
