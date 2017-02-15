% This generates the cylinders which are used to create the holes in the
% NEONATAL skull. The output is a FREECAD macro which creates ~3000 cylinders. This
% solid file is then removed from the skull model to create all the holes

%load the mesh created by Get_Approx_Mesh_For_Holes
load(['output' filesep 'Mesh_for_holes.mat']);

% find the points equally spaced on the surface, then find the normal
% vectors at those points. Then create a FreeCAD macro to generate them
generate_hole_locations_NN_skull(Mesh);