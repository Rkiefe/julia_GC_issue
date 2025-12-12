#=
    This is the code that runs a loop, where each iteration within
    that loop increases the memory consumption until an OOM issue

    !!! REQUIRES THE FOLLOWING PACKAGES !!!
    Gmsh, LinearAlgebra, SparseArrays, DelimitedFiles, Dierckx
=#

include("src/Femeko.jl")
include("src/magneticProperties.jl")

function MagnetoStaticSimulation(T::Float64, meshSize=0, localSize=0, showGmsh=true, verbose=false)
    #=
        Makes a large sparse matrix A and updates it inside a while loop
        Gmsh should make a mesh with about 300k elements and
        each while loop (Picard and Newton-Raphson) should take about 15 Gb
        to complete.

        Without GB.gc() inside the Picard loop
    =#

    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # Temperature
    # T::Float64 = 294.0

    # Applied field | T
    Hext::Vector{Float64} = [0.0, 0.0, 0.0]

    # Convergence criteria
    picardDeviation::Float64 = 1e-6
    maxDeviation::Float64 = Inf # Inf -> Don't run the N-R method
    maxAtt::Int32 = 500
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    data = DATA()

    # Load Gd
    density::Float64 = 7.9 # g/cm3
    loadMaterial( data,
                  "Materials",  # Folder with materials
                  "Gd_MFT",     # Data folder of target material
                  "Gd",         # Material name
                  density,      # g/cm3
                  T)

    # Spline for interpolation of the permeability
    spl = Spline1D(data.HofM, data.mu
                   ; bc="error" # nearest , extrapolate , error
                   )

    spl_dmu = Spline1D( data.HofM,
                        data.dmu
                      ; bc="error" # nearest , extrapolate , error
                      ) 

    # 3D Model
    gmsh.initialize()
    
    # Array of volume cell IDs
    cells = []

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells, true)

    # Add material
    addCuboid([0,0,0], [1.0, 1.2, 5.0], cells)
    box = addSphere([0,0,0], 25.0) # Add a bounding shell

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, false)
    
    # Print number of elements and nodes
    println("\nNumber of elements: ", mesh.nt)
    println("Number of internal elements: ", mesh.nInside)
    println("Number of internal nodes: ", mesh.nInsideNodes, "\n")

    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        centroids[:, k] = mean(mesh.p[:, nds], 2)
    end

    # Local stiffness matrix
    Ak::Matrix{Float64} = localStiffnessMatrix(mesh)

    # Convert to a compressed sparse column format
    rowIDs::Vector{Int} = zeros(16*mesh.nt)
    colIDs::Vector{Int} = zeros(16*mesh.nt)
    Acsc = zeros(length(rowIDs)) # Compressed sparse column format of the stiffness matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            for k in 1:mesh.nt
                rowIDs[(n-1)*mesh.nt + k] = mesh.t[i, k]
                colIDs[(n-1)*mesh.nt + k] = mesh.t[j, k]
            end
        end
    end # Compressed sparse stiffness matrix indices

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Magnetic permeability
    mu::Vector{Float64} = zeros(mesh.nt)

    # Hysteresis loop
    Hspan = vcat(range(0.0, 0.1, 5), range(0.12, 1.2, 10)) # , range(1.0, 0.0, 10)
    # Hspan = range(0.9, 1.0, 2)

    M_H::Vector{Float64} = zeros(length(Hspan)) 
    for (i, h) in enumerate(Hspan)
        println("\n", i, "/", length(Hspan))

        # Set the applied field
        Hext[3] = h # Applied field (T)

        # Boundary conditions
        RHS::Vector{Float64} = BoundaryIntegral(mesh, Hext, shell_id)

        # Reset the permeability
        mu .= mu0

        # FEM
        u::Vector{Float64} = zeros(mesh.nv+1)
        
        Hfield::Matrix{Float64} = zeros(3, mesh.nt)
        H::Vector{Float64} = zeros(mesh.nt)
        Hold::Vector{Float64} = zeros(mesh.nt)

        att::Int32 = 0
        div::Float64 = Inf
        while div > picardDeviation && att < maxAtt 

            att += 1
            Hold .= H

            # Updated the compressed stiffness matrix
            for n in 1:16
                for k in 1:mesh.nt
                    Acsc[(n-1)*mesh.nt + k] = Ak[n, k] * mu[k]
                end
            end # Acsc

            # Update the global stiffness matrix
            A = sparse(rowIDs, colIDs, Acsc, mesh.nv, mesh.nv)

            # Magnetic scalar potential
            u = [A Lag;Lag' 0]\[-RHS;0]
            
            # Check solution
            if any(x -> !isfinite(x), u)
                error("Nans/Infs in the scalar potential")
            end

            # Magnetic field
            Hfield .= 0.0
            for k in 1:mesh.nt
                nds = mesh.t[:,k];

                # Sum the contributions
                for nd in nds
                    # obtain the element parameters
                    _, b, c, d = abcd(mesh.p,nds,nd)

                    Hfield[1, k] -= u[nd]*b;
                    Hfield[2, k] -= u[nd]*c;
                    Hfield[3, k] -= u[nd]*d;
                end
            end

            # Magnetic field intensity
            for k in 1:mesh.nt
                H[k] = norm(Hfield[:, k])
            end

            # Update magnetic permeability            
            mu[mesh.InsideElements] .= spl(H[mesh.InsideElements])

            # Check deviation from previous result
            div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
            verbose ? println(att, " | mu0 |H(n)-H(n-1)| = ", div) : nothing

            # Ask for garbage collection because otherwise there is a memory build up
            # Still don't know why
            # GC.gc()

        end # Picard Iteration

        # Magnetization
        chi::Vector{Float64} = mu./mu0 .- 1.0

        Mfield::Matrix{Float64} = zeros(3, mesh.nt)
        Mfield[1, :] = chi.*Hfield[1, :]
        Mfield[2, :] = chi.*Hfield[2, :]
        Mfield[3, :] = chi.*Hfield[3, :]
        
        # Average magnetization
        M_avg::Float64 = 0.0
        volume::Float64 = 0.0
        for k in mesh.InsideElements
            volume += mesh.VE[k]
            M_avg  += mesh.VE[k]*Mfield[3, k]
        end # <H>

        # Store the M(H)
        M_H[i] = M_avg/volume
    end

end # end of main

meshSize  = 2.0
localSize = 0.05
showGmsh = false
verbose = true

# For each temperature, run a simulation
for T in 270.0:4.0:330.0
    MagnetoStaticSimulation(T, meshSize, localSize, showGmsh, verbose)
end
