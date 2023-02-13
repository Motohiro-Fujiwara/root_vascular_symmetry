# Root vascular symmetry
 Boundary quantification and vertex model code for root vascular symmetry in "Patterned proliferation orients tissue-wide stress to control root vascular symmetry in Arabidopsis" (https://doi.org/10.1016/j.cub.2023.01.036). Description of the contents of each folder.
 1.	boundary_quantification
<br>We determined each position of the tricellular junctions in the root vascular tissue at the transverse section of the proliferation zone using the Fiji plugin Tissue Analyzer.
<br>Symmetry is calculated by boundary_symmetry.xlsx.
<br>Roughness is calculated by boundary_roughness.xlsx.
<br>Angle is calculated by boundary_angle.xlsx.
<br>Aspect ratio is calculated by Fiji plugin Tissue Analyzer.
2.	vertex_model
<br>Initial cell position data set is written in initial_cell_position folder.
<br>Calculation 
<br>C++:cellform_vascul_growth.cpp
<br>Visualize
<br>Processing:sketch_vascul_sim.pde
