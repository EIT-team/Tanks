% requires iso2mesh and stlTools https://uk.mathworks.com/matlabcentral/fileexchange/51200-stltools


%%
% resolution of output volumetric image
vol_res=0.5; % size of voxels in mm
pixel_scale = 1/vol_res;


%% Loading stls
disp('Loading Stl meshes');
[skull.vertices,skull.faces,skull.normals,skull.name] = stlRead('Bodies_For_Meshing_Skull.stl');
[scalp.vertices,scalp.faces,scalp.normals,scalp.name] = stlRead('Bodies_For_Meshing_Scalp.stl');
stlPlot(skull.vertices,skull.faces,skull.name);
stlPlot(scalp.vertices,scalp.faces,scalp.name);

figure;
hold on
plotmesh(scalp.vertices,scalp.faces)
plotmesh(skull.vertices,skull.faces)
hold off

drawnow
%% convert to binary mask

disp('Making binary masks')

%define the mask range - this forces all the masks to be the same size
mask_range= [floor(min([skull.vertices; scalp.vertices]))-vol_res; ceil(max([skull.vertices; scalp.vertices]))+vol_res];

%find the middle for plotting the slices
mask_mid_idx = round(((mask_range(2,:) - mask_range(1,:))/2)/vol_res) ;
mask_mid = mask_mid_idx + mask_range(1,:);

%convert the skull mask
[skull.mask, skull.transform]=surf2vol(skull.vertices,skull.faces,mask_range(1,1):vol_res:mask_range(2,1),mask_range(1,2):vol_res:mask_range(2,2),mask_range(1,3):vol_res:mask_range(2,3),'fill',1);

%convert scalp mask
[scalp.mask, scalp.transform] =surf2vol(scalp.vertices,scalp.faces,mask_range(1,1):vol_res:mask_range(2,1),mask_range(1,2):vol_res:mask_range(2,2),mask_range(1,3):vol_res:mask_range(2,3),'fill',1);

figure
imagesc(skull.mask(:,:,mask_mid_idx(3)));daspect([1,1,1])
title('Skull mask')

figure
imagesc(scalp.mask(:,:,mask_mid_idx(3)));daspect([1,1,1])
title('Scalp mask');
% combine masks
full_mask = scalp.mask + skull.mask;

figure
imagesc(full_mask(:,:,mask_mid_idx(3)));daspect([1,1,1])
title('Combined mask')

drawnow

%% convert positions
elec_pos=dlmread('NNelecposorig.txt');

figure
hold on
plotmesh(scalp.vertices,scalp.faces)
plot3(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3),'o');
hold off
title('Elecs in original positions');

%adjust positons based on the affine transform output from surf2vol
newx=elec_pos(:,1) -scalp.transform(1,4);
newy=elec_pos(:,2) -scalp.transform(2,4);
newz=elec_pos(:,3) -scalp.transform(3,4);

elec_pos_new=[newx newy newz];

%% save stuff for mesher

%save the volumetric data
saveinr_EIT(uint8(full_mask),'NNvol.inr',vol_res*[1 1 1]);
saveinr_EIT(uint8(scalp.mask),'NNvol_homo.inr',vol_res*[1 1 1]);

% save the electrode locations
dlmwrite('NNelecpos.txt',elec_pos_new*pixel_scale);

%% RUN MESHER

%% load mesher output

% Mesh=load('Mesh.mat');
% figure
% hold on
% DisplayBoundaries(Mesh)
% % plot3(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3),'o');
% plot3(elec_pos_new(:,1),elec_pos_new(:,2),elec_pos_new(:,3),'r.','Markersize',30);
% hold off







