% This converts the comsol mesh generated from the Solid Model, into the
% "Mesh" structure for use in EIT software. Requires COMSOL LIVE LINK
%
% This is for the complete Adult mesh for the forward solver

model=mphload('AdultTank_Mesh_4mln.mph');
[stats,Data]=mphmeshstats(model);

tri=Data.elem{2};
tri=double(tri');
Mesh.Tetra=double(tri)+ones(size(tri,1),size(tri,2));
Mesh.Nodes=Data.vertex';
Mesh.mat_ref=double(Data.elementity{2});

%save
save(['output' filesep 'AdultTank_Mesh_4mln.mat'],'Mesh','-v7.3');