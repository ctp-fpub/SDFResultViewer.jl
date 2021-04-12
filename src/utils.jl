function sdfcount(root_folder, pattern)
    folders = filter!(isdir, readdir(root_folder, join=true))
    filter!(f-> occursin(pattern, f), folders)
    n = [count(f->endswith(f, ".sdf"), readdir(dir)) for dir in folders]
    idx = findall(!iszero, n)

    Plots.plot(basename.(folders[idx]), n[idx], rotation=45, legend=false)
end

function similar_laser(laser::Type{LaguerreGaussLaser}, t_profile, file)
    fx = get_parameter(file, :constant, :f_x)*u"m" |> unit_l
    λ = get_parameter(file, :laser, :lambda)
    m = get_parameter(file, :constant, :m)
    a₀ = float(get_parameter(file, :constant, :a0))
    w₀ = get_parameter(file, :constant, :w_0)u"m" |> unit_l
    τ = get_parameter(file, :constant, :tau_0)u"m" |> unit_t
    t = get_time(file) |> unit_t

    s = setup_laser(laser, :SI_unitful;
        λ,
        m,
        p = 0,
        a₀,
        w₀,
        # ϕ₀=π/2,   # aparently this doesn't do anything
        profile = t_profile(;
            c = c_0,
            τ,
            z₀ = -fx,
            t₀ = fx/(π*c_0) |> unit_t
        )
    )
end

function field_rotation(dir)
    if dir == :x
        RotZ(π)*RotY(-π/2)*RotZ(-π/2)
    elseif dir == :y
        RotY(-π/2)
    else
        1
    end
end

function similar_E(::Val{true}, laser, t, propagation_dir, component::Int, grid; x₀, y₀, z₀)
    map(Iterators.product(grid...)) do (x,y,z)
        r = SVector{3}(x-x₀, y-y₀, z-z₀)
        R = field_rotation(propagation_dir)
        analytic_E = R * (LaserTypes.E(R \ r, laser) * LaserTypes.g((R \ r)[3], t, laser))
        analytic_E[component]
    end
end

function similar_E(::Val{false}, laser, t, propagation_dir, component::Int, grid; x₀, y₀, z₀)
    map(Iterators.product(grid...)) do (x,y,z)
        r = SVector{3}(x-x₀, y-y₀, z-z₀)
        R = field_rotation(propagation_dir)
        analytic_E = R * LaserTypes.E(R \ r, t, laser)
        analytic_E[component]
    end
end

function similar_E(E, t, laser;
    component=1,
    propagation_dir=:x,
    complex=false,
    x₀=zero(recursive_bottom_eltype(getdomain(E))),
    y₀=zero(recursive_bottom_eltype(getdomain(E))),
    z₀=zero(recursive_bottom_eltype(getdomain(E))))

    E_approx = approximate_field(E)

    grid = AxisGrid((getdomain(E_approx).*unit_l...,))
    data = similar_E(Val(complex), laser, t, propagation_dir, component, grid; x₀,y₀,z₀)

    ScalarField(data, grid)
end

function compare_E_slice_with_analytic(laser, profile, file, dir;
    slice_location=0unit_l,
    slice_dir=:x)

    fx = get_parameter(file, :constant, :f_x)*u"m" |> unit_l
    t = get_time(file) |> unit_t
    s = similar_laser(laser, profile, file)

    E = uconvert(unit_E, unit_l, file["e$dir"])

    analytic_E = similar_E(E, t, s,
        component=dir_to_idx(dir),
        propagation_dir=:x,
        x₀=fx) |> unit_E

    plt1 = E_slice_plot(analytic_E, dir, slice_location, t;
        legend_title = "Analytic",
        slice_dir,
        issliced=false)
    plt2 = E_slice_plot(E, dir, slice_location, t;
        legend_title = "Numeric",
        slice_dir,
        issliced=false)
    Plots.plot(plt1, plt2, size=(400,600), layout=(2,1))
end
