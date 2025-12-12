# julia_GC_issue
A repo to replicate an OOM issue with Julia

## Must add the followig Julia packages!!
```julia
add Gmsh, LinearAlgebra, SparseArrays, DelimitedFiles, Dierckx
```

Run the main.jl and it will lead to an OOM. If you uncomment the `GC.gc()` on line 201 it will not.
