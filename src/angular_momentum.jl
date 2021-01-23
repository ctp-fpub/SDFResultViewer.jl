function compute_L(file, species)
    w, (x,y,z), px, py, pz = file["weight/"*species,
                                  "grid/"*species,
                                  "px/"*species,
                                  "py/"*species,
                                  "pz/"*species]

    r = [SVector{3}(x[i],y[i],z[i]) for i in eachindex(x)]
    p = VectorVariable(px, py, pz)

    L = VectorOfArray(w .* r .× p)
end

function compute_Lx(file, species)
    w, r, py, pz = file["weight/"*species,
                        "grid/"*species,
                        "py/"*species,
                        "pz/"*species]

    y = r[2]
    z = r[3]
    Lx = @. w * (y * pz - z * py)
end
