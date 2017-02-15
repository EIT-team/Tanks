% This generates the cylinders which are used to create the holes in the
% skull. The output is a FREECAD macro which creates ~3000 cylinders. This
% solid file is then removed from the skull model to create all the holes

%load the mesh created by Get_Approx_Mesh_For_Holes
load(['Mesh_for_construction.mat']);

% find the points equally spaced on the surface, then find the normal
% vectors at those points. adjust the diameter based on the desired
% distribution. Then create a FreeCAD macro to generate them
generate_holes_locations_Adult_skull(Mesh);