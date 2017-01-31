% This converts the comsol mesh generated from the Solid Model, into the
% "Mesh" structure for use in EIT software. Requires COMSOL LIVE LINK
%
% This is the approximation of the surface used to generate the surface
% normal vectors used in generating the cylinders in the skull

% load the approximate SURFACE mesh
model=mphload('for_holes.mph');
[stats,Data]=mphmeshstats(model);

% load the surface triangle
Tri=Data.elem{2};
Tri=double(Tri');

%Generate Holes function wants a little different values
Mesh.Faces=double(Tri)+ones(size(Tri,1),size(Tri,2));
Mesh.Node_faces=1000*Data.vertex';
Mesh.mat_ref=double(Data.elementity{2});
save(['output' filesep 'Mesh_for_holes.mat'],'Mesh','-v7.3');