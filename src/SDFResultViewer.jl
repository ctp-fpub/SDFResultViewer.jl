module SDFResultViewer

export compute_L, compute_Lx, mean_Lx, Lx_sclice, subsample,
    read_positions, read_E, get_times

using SDFReader
using RecursiveArrayTools
using LinearAlgebra
using Unitful
using PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
using StaticArrays
using ImageTransformations: imresize
using IntervalSets
using FileTrees
using Transducers, ThreadsX

function compute_L(file)
    keys = ["weight/electron",
            "grid/electron",
            "px/electron",
            "py/electron",
            "pz/electron"]

    w, (x,y,z), px, py, pz = readkeys(file, keys)

    r = [SVector{3}(x[i],y[i],z[i]) for i in eachindex(x)]
    p = [SVector{3}(px[i],py[i],pz[i]) for i in eachindex(px)]

    L = VectorOfArray(w .* r .× p)
end

function compute_Lx(file)
    keys = ["weight/electron",
            "grid/electron",
            "py/electron",
            "pz/electron"]

    w, (x,y,z), py, pz = readkeys(file, keys)

    Lx = @. w * (y * pz - z * py)
end

function compute_Lx(file, λ)
    ω = 2π*c_0/λ
    unit_L = m_e * c_0^2 / ω

    Lx = uconvert.(NoUnits, compute_Lx(file) / unit_L)
end

function read_E(file)
    k = ["ex", "ey", "ez"]
    blocks = file_summary(file)

    grid, Ex, Ey, Ez = readkeys(file, k)

    E = SVector{3}.(Ex, Ey, Ez)

    return grid, E
end

function subsample(v, target_size)
    length(v) > target_size ? imresize(v, target_size) : v
end

function read_positions(file, λ)
    ω = 2π*c_0/λ
    unit_length = c_0 / ω
    keys = ["grid/electron"]

    ((x,y,z),) = readkeys(file, keys)

    return uconvert.(NoUnits, x / unit_length),
           uconvert.(NoUnits, y / unit_length),
           uconvert.(NoUnits, z / unit_length)
end

function sdf_tree(dir)
    paths = filter(f->endswith(f, ".sdf"), readdir(dir))
    maketree(dir=>paths)
end

function mean_Lx_dag(dir, λ)
    tree = sdf_tree(dir)

    Lx = FileTrees.load(tree, lazy=true) do file
        println("Loading $(path(file)) on thread $(Threads.threadid())")
        compute_Lx(path(file), λ)
    end
    exec(mapvalues(mean, Lx))
end

function mean_Lx(dir, λ; cond=lx->true)
    paths = filter(f->endswith(f, ".sdf"), readdir(dir))

    ThreadsX.map(paths) do f
        file = joinpath(dir, f)
        println("Loading $file")
        Lx = compute_Lx(file, λ)
        z = zero(eltype(Lx))
        (Lx ./ length(Lx)) |>
            Filter(cond) |>
            foldxt(+, simd=true, init=z)
    end
end

function Lx_sclice(file, λ, slice_location, ϵ=2e-3; target_size=3*10^4)
    x,y,z = read_positions(file, λ)
    ids = filter(i-> x[i] ∈ slice_location ± ϵ, axes(x, 1))

    Lx = subsample_Lx(file, λ, ids, target_size)

    Lx, subsample(y[ids], target_size), subsample(z[ids], target_size)
end

function subsample_Lx(file, λ, id, target_size=3*10^4)
    Lx = compute_Lx(file, λ)
    subsample(Lx[id], target_size)
end

function get_times(dir, λ)
    paths = filter(f->endswith(f, ".sdf"), readdir(dir))
    ts = [get_time(joinpath(dir, f)) for f in paths]
    unit_t = λ / (2π*c_0)

    uconvert.(NoUnits, ts / unit_t)
end

end
