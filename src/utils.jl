function sdfcount(root_folder, pattern)
    folders = filter!(isdir, readdir(root_folder, join=true))
    filter!(f-> occursin(pattern, f), folders)
    n = [count(f->endswith(f, ".sdf"), readdir(dir)) for dir in folders]
    idx = findall(!iszero, n)

    Plots.plot(basename.(folders[idx]), n[idx], rotation=45, legend=false)
end
