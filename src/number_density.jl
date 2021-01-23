function average_linear_density_x(file, species)
    n = file["number_density/"*species]

    ω = 2π*c_0/λ
    unit_length = c_0 / ω

    x = uconvert.(NoUnits, n.grid[1]./unit_length)
    y = uconvert.(NoUnits, n.grid[2] ./ unit_length)
    z = uconvert.(NoUnits, n.grid[3] ./ unit_length)
    n_normed = uconvert.(NoUnits, n.data .* unit_length^3)
    nlx = ThreadsX.map(axes(n.data, 1)) do i
        integrate((y, z), view(n_normed, i,:,:))
    end

    x, nlx
end

function average_radial_density(file, λ)
    blocks = file_summary(file)

    n = open(file) do f
        SDFReader.read_scalar_field(f, blocks, "number_density/electron")
    end

    ω = 2π*c_0/λ
    unit_length = c_0 / ω

    x = uconvert.(NoUnits, n.grid[1]./unit_length)
    y = uconvert.(NoUnits, n.grid[2] ./ unit_length)
    z = uconvert.(NoUnits, n.grid[3] ./ unit_length)
    r = SVector.(Iterators.product(x,y,z))
    n_normed = uconvert.(NoUnits, n.data .* unit_length^3)

    coords = Cylindrical.(SVector{3}.(y,z,x))
    # n_cyl =

    nlx = ThreadsX.map(axes(n.data, 1)) do i
        integrate((y, z), view(n_normed, i,:,:))
    end

    x, nlx
end
