% This converts the comsol mesh generated from the Solid Model, into the
% "Mesh" structure for use in EIT software. Requires COMSOL LIVE LINK
%
% This is the approximation of the surface used to generate the surface
% normal vectors used in generating the cylinders in the skull

%load the mesh created by Get_Approx_Mesh_For_Holes
load(['output' filesep 'Mesh_for_holes.mat']);

% find the points equally spaced on the surface, then find the normal
% vectors at those points. Then create a FreeCAD macro to generate them
generate_hole_locations_NN_skull(Mesh);