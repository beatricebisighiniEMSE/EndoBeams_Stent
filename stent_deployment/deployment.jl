function deployment(initial_positions_stent_mat, initial_positions_stent, connectivity_stent, filename_sdf, output_dir_crimping, output_dir_positioning, output_dir_deployment)

    if !isdir(output_dir_deployment)
        println("Create directory")
        mkdir(output_dir_deployment)
    end

    dir = pwd()
    cd(output_dir_deployment)
    foreach(rm, filter(endswith(".txt"), readdir()))
    foreach(rm, filter(endswith(".vtk"), readdir())) 
    cd(dir)

    # -----------------------------------------------------------------------------------------
    # Stent 
    # -------------------------------------------------------------------------------------------

    positioned_positions_stent = readdlm(output_dir_positioning * "crimped.txt") + readdlm(output_dir_positioning * "u.txt")

    # -----------------------------------------------------------------------------------------
    # Building the nodes
    # -------------------------------------------------------------------------------------------

    # total number of nodes
    nnodes = size(initial_positions_stent_mat, 1)

    # initial conditions
    u⁰ = positioned_positions_stent - initial_positions_stent_mat
    u̇⁰ = zeros(nnodes, 3)
    ü⁰ = zeros(nnodes, 3)
    w⁰ = zeros(nnodes, 3)
    ẇ⁰ = zeros(nnodes, 3)
    ẅ⁰ = zeros(nnodes, 3)

    R⁰ = readdlm(output_dir_positioning * "R.txt") 
 
    # nodes StructArray
    nodes = build_nodes(initial_positions_stent_mat, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, R⁰)

    # -------------------------------------------------------------------------------------------
    # Building the beams
    # -------------------------------------------------------------------------------------------

    # geometric and material properties
    E = 206*1e3
    ν = 0.33
    ρ = 8*1e-9 
    damping = 1E4*5

    Re₀ = read_ics_mat(readdlm(output_dir_crimping * "Re0.txt"))

    # beams vector
    beams = build_beams(nodes, connectivity_stent, E, ν, ρ, rWireSection, damping, Re₀)

    # contact parameters
    kₙ = 4/3 * 5/(1-0.5^2)*sqrt(rWireSection) # Approximate Hertz contact with 5 MPa wall stiffness
    μ = 0.5
    εᵗ = 0.1 #regularized parameter for friction contact
    ηₙ = 0.1
    kₜ = kₙ
    ηₜ = ηₙ
    u̇ₛ = 0.005
    contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)

    # -------------------------------------------------------------------------------------------
    # External forces
    # -------------------------------------------------------------------------------------------

    # external force and applied dof
    loaded_dofs = []
    force(t, node_idx) = 0

    ext_forces = ExternalForces(force, loaded_dofs)

    # -------------------------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------------------------

    # number of dof (6 per node)
    ndofs = nnodes*6

    # penalty constraints
    kᶜᵒⁿ = kₙ
    ηᶜᵒⁿ = kₙ
    nodespairs = get_nodespairs_stent(initial_positions_stent)
    writedlm("nodespairs.txt", nodespairs)
    constraints = build_constraints(nodespairs, kᶜᵒⁿ, ηᶜᵒⁿ)

    # Dirichlet boundary conditions: blocked positions
    fixed_dofs = Int[]
    free_dofs = setdiff(1:ndofs, fixed_dofs)

    # Dirichlet dof (x6)
    disp_dofs = Int[]
    disp_vals = Float64[]
    disp(t, node_idx) = 0

    # boundary conditions strucutre
    bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)

    # -------------------------------------------------------------------------------------------
    # SDF
    # -------------------------------------------------------------------------------------------

    sdf = Discrete_SDF(filename_sdf, rWireSection, true)

    # -------------------------------------------------------------------------------------------
    # Final configuration
    # -------------------------------------------------------------------------------------------

    # configuration: mesh, external forces and boundary conditions
    conf = Configuration(nodes, beams, constraints, ext_forces, bcs, contact, sdf)

    # -------------------------------------------------------------------------------------------
    # Time stepping parameters
    # -------------------------------------------------------------------------------------------

    # initial time step and total time
    ini_Δt = 1e-9
    max_Δt = 1e-5
    Δt_plot = 1e-5
    tᵉⁿᵈ = 100

    params = Params(ini_Δt = ini_Δt, Δt_plot = Δt_plot, max_Δt = max_Δt, tᵉⁿᵈ = tᵉⁿᵈ, output_dir = output_dir_deployment, stop_on_energy_threshold=true, energy_threshold=5*1e-5, tol_res = 1e-3, tol_ΔD = 1e-3, stop_long_simulation=false, record_timings=false, verbose=true)

    # -------------------------------------------------------------------------------------------z
    # Start simulation
    # -------------------------------------------------------------------------------------------

    solver!(conf, params)
    write_txt_solution(nodes, beams, length(nodes), length(beams), output_dir_deployment, false)

end 

