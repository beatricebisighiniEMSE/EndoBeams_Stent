using DelimitedFiles
using Parameters
using LinearAlgebra
using Dierckx
using StaticArrays
using LinearAlgebra
using StructArrays
using EndoBeams

include("utils_stent.jl")
include("crimping.jl")
include("positioning_cl.jl")
include("positioning.jl")
include("deployment.jl")

# -------------------------------------------------------------------------------------------
# Create stent 
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStent()
initial_positions_stent, connectivity_stent = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
initial_positions_stent_mat = reshape(reinterpret(Float64, initial_positions_stent), (3, length(initial_positions_stent)))'

writedlm("pos.txt", initial_positions_stent)
writedlm("conn.txt", connectivity_stent)

# -------------------------------------------------------------------------------------------
# Select files 
# -------------------------------------------------------------------------------------------

i = 0 # which case 
filename_surf = "stent_deployment/input/mesh_$i.stl"
filename_cl = "stent_deployment/input/cl_$i.vtk"
filename_sdf = "stent_deployment/input/sdf_$i.vtk"

# -------------------------------------------------------------------------------------------
# Crimping
# -------------------------------------------------------------------------------------------

output_dir_crimping = "stent_deployment/output3D/outputCrimping3D/"
crimping(rStent, rCrimpedStent, initial_positions_stent, connectivity_stent, output_dir_crimping)

output_dir_positioning_cl = "stent_deployment/output3D/outputPositioningCl3D/"
output_dir_positioning = "stent_deployment/output3D/outputPositioning3D/"
output_dir_deployment = "stent_deployment/output3D/outputDeployment3D/"

# -------------------------------------------------------------------------------------------
# Geometrical positioning centerline
# -------------------------------------------------------------------------------------------


nb_iterations = 20 #how many interations for the positioning
deploy_pos = get_init_pos_deploy_middle(filename_cl, 0.7, initial_positions_stent, output_dir_crimping)

initial_cl_stent = positioning_cl(initial_positions_stent, connectivity_stent, nb_iterations, deploy_pos, filename_cl, output_dir_crimping, output_dir_positioning_cl)

# -------------------------------------------------------------------------------------------
# Physical positioning stent_deployment
# -------------------------------------------------------------------------------------------

positioning(initial_positions_stent, connectivity_stent, nb_iterations, initial_cl_stent[1], output_dir_crimping, output_dir_positioning_cl, output_dir_positioning)

crimped_positions_stent = initial_positions_stent .+  read_ics_vec(readdlm(output_dir_crimping * "u.txt"))
positions_cl =  get_centerline_stent(crimped_positions_stent)
set_origin!(crimped_positions_stent, positions_cl[1]-initial_cl_stent[1])
open(output_dir_positioning * "/crimped.txt", "w") do io
    writedlm(io, crimped_positions_stent)
end

disp = read_ics_vec(readdlm(output_dir_positioning * "u.txt"))
write_vtk_configuration("stent_deployment/output3D/" * "positioned_$i.vtk", crimped_positions_stent + disp, connectivity_stent)

# -------------------------------------------------------------------------------------------
# Deployment
# -------------------------------------------------------------------------------------------

deployment(initial_positions_stent_mat, initial_positions_stent, connectivity_stent, filename_sdf, output_dir_crimping, output_dir_positioning, output_dir_deployment)

disp = read_ics_vec(readdlm(output_dir_deployment * "u.txt"))
write_vtk_configuration("stent_deployment/output3D/" * "deployed_$i.vtk", initial_positions_stent + disp, connectivity_stent)
