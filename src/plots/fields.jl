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

function compare_slice_with_analytic(numeric_f, f, laser, profile, file;
    component=1,
    slice_location=0unit_l,
    slice_dir=:x,
    label="",
    f_unit)

    fx = get_parameter(file, :constant, :f_x)*u"m" |> unit_l
    t = get_time(file) |> unit_t
    s = analytic_laser(laser, profile, file)

    grid = getdomain(numeric_f)

    analytic_f = apply_analytic(f, grid, t, s;
        component,
        propagation_dir=:x,
        x₀=fx) |> f_unit

    fig = Figure(resolution = (800, 700))

    units = " (" * string(unit(recursive_bottom_eltype(grid))) * ")"
    xlabel, ylabel = filter(x->x≠string(slice_dir),
        ["x", "y", "z"] .* units)
    # ts = @lift sprint(show, $t, context=:compact => true)
    # loc = @lift sprint(show, $slice_location, context=:compact => true)

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

    analytic_slice, loc = to_observable(analytic_f; slice_dir, slice_location, issliced=false)
    numeric_slice, loc = to_observable(numeric_f; slice_dir, slice_location, issliced=false)

    plt1 = fieldplot!(ax1, analytic_slice)
    plt2 = fieldplot!(ax2, numeric_slice)

    Colorbar(fig[1,1][2,1], plt1;
        width = Relative(3/4),
        height = 10,
        vertical = false,
        tellheight = true,
        label,
    )
    Colorbar(fig[1,2][2,1], plt2;
        width = Relative(3/4),
        height = 10,
        vertical = false,
        tellheight = true,
        label,
    )

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

    return plt
end
