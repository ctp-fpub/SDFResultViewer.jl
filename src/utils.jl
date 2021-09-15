function sdfcount(root_folder, pattern)
    folders = filter!(isdir, readdir(root_folder, join=true))
    filter!(f-> occursin(pattern, f), folders)
    n = [count(f->endswith(f, ".sdf"), readdir(dir)) for dir in folders]
    idx = findall(!iszero, n)

    Plots.plot(basename.(folders[idx]), n[idx], rotation=45, legend=false)
end

function read_EM_fields(file)
    Ex, Ey, Ez = uconvert.(unit_E, unit_l, file[:ex, :ey, :ez])
    Bx, By, Bz = uconvert.(unit_B, unit_l, file[:bx, :by, :bz])

    E = build_vector((Ex, Ey, Ez), (:x, :y, :z))
    B = build_vector((Bx, By, Bz), (:x, :y, :z))

    E, B
end

function dynamic_slice(f, slice_dir, slice_location)
    if f isa Observable
        obs_f = ustrip(f)
    else
        @debug "Converted field to Observable"
        obs_f = Observable(ustrip(f))
    end
    if slice_location isa Observable
        loc = @lift ustrip($slice_location)
    else
        @debug "Converted slice_location to Observable"
        loc = Observable(ustrip(slice_location))
    end

    f_slice = @lift slice($obs_f, slice_dir, $loc)

    return f_slice, loc
end

function analytic_laser(laser::Type{LaguerreGaussLaser}, t_profile, file)
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

function apply_analytic(::Val{true}, f, laser, t, propagation_dir, component, grid; x₀, y₀, z₀)
    mapgrid(grid) do (x,y,z)
        r = SVector{3}(x-x₀, y-y₀, z-z₀)
        R = field_rotation(propagation_dir)
        analytic_E = R * (f(R \ r, laser) * LaserTypes.g((R \ r)[3], t, laser))

        if component isa Number
            analytic_E[component]
        else
            analytic_E
        end
    end
end

function apply_analytic(::Val{false}, f, laser, t, propagation_dir, component, grid; x₀, y₀, z₀)
    mapgrid(grid) do (x,y,z)
        r = SVector{3}(x-x₀, y-y₀, z-z₀)
        R = field_rotation(propagation_dir)
        analytic_E = R * f(R \ r, t, laser)

        if component isa Number
            analytic_E[component]
        else
            analytic_E
        end
    end
end

function apply_analytic(f, grid, t, laser;
    component=1,
    propagation_dir=:x,
    complex=false,
    downsample_grid=true,
    x₀=zero(recursive_bottom_eltype(grid)),
    y₀=zero(recursive_bottom_eltype(grid)),
    z₀=zero(recursive_bottom_eltype(grid)))

    grid = downsample_grid ? downsample_approx(grid) : grid
    data = apply_analytic(Val(complex),
        f,
        laser,
        t,
        propagation_dir,
        component,
        grid;
        x₀,y₀,z₀)
    if eltype(data) <: Number
        ScalarField(data, grid)
    else
        # broken
        @error "Use analytic_E or analytic_B"
        VectorField(data, grid)
    end
end

# Stopgap solution untill we can greate vector fields correctly with the above
function get_analytic_E(grid, t, s;
    downsample_grid=true,
    x₀=zero(recursive_bottom_eltype(grid)),
    y₀=zero(recursive_bottom_eltype(grid)),
    z₀=zero(recursive_bottom_eltype(grid)))

    grid = downsample_grid ? downsample_approx(grid) : grid

    analytic_Ex = apply_analytic(LaserTypes.E, grid, t, s;
        x₀,y₀,z₀, component=1) |> unit_E
    analytic_Ey = apply_analytic(LaserTypes.E, grid, t, s;
        x₀,y₀,z₀, component=2) |> unit_E
    analytic_Ez = apply_analytic(LaserTypes.E, grid, t, s;
        x₀,y₀,z₀, component=3) |> unit_E

    build_vector((analytic_Ex, analytic_Ey, analytic_Ez), (:x,:y,:z))
end

function get_analytic_B(grid, t, s;
    downsample_grid=true,
    x₀=zero(recursive_bottom_eltype(grid)),
    y₀=zero(recursive_bottom_eltype(grid)),
    z₀=zero(recursive_bottom_eltype(grid)))

    grid = downsample_grid ? downsample_approx(grid) : grid

    analytic_Bx = apply_analytic(LaserTypes.B, grid, t, s;
        x₀,y₀,z₀, component=1) |> unit_B
    analytic_By = apply_analytic(LaserTypes.B, grid, t, s;
        x₀,y₀,z₀, component=2) |> unit_B
    analytic_Bz = apply_analytic(LaserTypes.B, grid, t, s;
        x₀,y₀,z₀, component=3) |> unit_B

    build_vector((analytic_Bx, analytic_By, analytic_Bz), (:x,:y,:z))
end
