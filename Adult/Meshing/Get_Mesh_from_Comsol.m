% This converts the comsol mesh generated from the Solid Model, into the
% "Mesh" structure for use in EIT software. Requires COMSOL LIVE LINK
%
% This is for the complete Adult mesh used in reconstructions

%load the tetrahedral Neonatal head mesh
model=mphload('AdultTank_for_construction.mph');
[meshstats,meshdata] = mphmeshstats(model,'mesh1');

%load the tetrahedra and rearrange into correct format
Mesh.Faces=(double(meshdata.elem{3})+ones(size(meshdata.elem{3})))';
Mesh.Nodes=meshdata.vertex';
Mesh.Tetra=(double(meshdata.elem{2})+ones(size(meshdata.elem{2})))';
Mesh.Nodes=Mesh.Nodes(:,[1,3,2])*1000;
% remove any isolated nodes
[Mesh.Nodes_faces, Mesh.Faces]=removeisolatednode_tri(Mesh.Nodes(:,1:3),Mesh.Faces);
Mesh.mat_ref=meshdata.elementity{2};
%save
save(['output' filesep 'Mesh_for_construction.mat'],'Mesh','-v7.3');