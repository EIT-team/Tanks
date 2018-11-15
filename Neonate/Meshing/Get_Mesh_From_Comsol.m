% This converts the comsol mesh generated from the Solid Model, into the
% "Mesh" structure for use in EIT software. Requires COMSOL LIVE LINK
%
% This is for the complete Neonatal mesh used in reconstructions
%
% THIS IS NOW OBSOLETE, USE SurfToSegmentation INSTEAD

%load the tetrahedral Neonatal head mesh
model=mphload('4p5mln_mesh.mph');
[stats,Data]=mphmeshstats(model);

%load the tetrahedra
Tetra=Data.elem{2};
Tetra=double(Tetra');
%rearrange into correct format
Mesh.Tetra=double(Tetra)+ones(size(Tetra,1),size(Tetra,2));
Mesh.Nodes=1000*Data.vertex';
Mesh.mat_ref=double(Data.elementity{2});
%save
save(['output' filesep 'Mesh_4p5_mln_forward.mat'],'Mesh','-v7.3');
