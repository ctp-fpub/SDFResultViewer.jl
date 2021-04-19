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
