
%filepath=[filepath '/'];
filepath=[pwd filesep 'output' filesep];

disp('Loading Mesh');
%load mesh
load('AdultTank_Mesh_4mln.mat');
%load electrode positions
load('Adult_Tank_Pos.mat');

%% Ground node
disp('Finding Ground node');
%find ground node as node at bottom
gnd_pos=Mesh.Nodes(Mesh.Nodes(:,3)==min(Mesh.Nodes(:,3)),:);
gnd_pos=gnd_pos(1,1:3);

%% Protocol

disp('Finding Protocol');
%find "optimal" protocol
Prot_Adult16=dlmread('NCopt_injection_protocol.txt');


%find complete prt
[Adult_Prt_full,Adult_keep_idx,Adult_rem_idx]=ScouseTom_data_findprt(Prot_Adult16,32); %from ScouseTom Repo

save('Adult2016Prt','Adult_keep_idx','Adult_Prt_full'); %for ref
save([filepath 'Adult2016Prt'],'Adult_keep_idx','Adult_Prt_full');

%write protocol to textfile for forward solver
dlmwrite([filepath 'Adult2016_Prt_full.txt'],Adult_Prt_full);

%write protocol injections to textfile for ScouseTom System
dlmwrite([filepath 'Adult2016_Prt_Injs.txt'],Prot_Adult16);

%% ExpSetup

disp('Sorting ExpSetup');
%load previously used with with old protocol
load(['Example_ExpSetup'],'ExpSetup');

%first make ExpSetup for Tank Comp Settings - this is at 300uA at 1khz. We
%cannot use this to make images, but we can use it to check tank is ok

%change protocol to new one
ExpSetup.Protocol=Prot_Adult16;

ExpSetup.Elec_num=32;
ExpSetup.Amp=240;
ExpSetup.Freq=1725;

%validate the new settings - ignore warning regarding too high current
[okflag,ExpSetup]=ScouseTom_ValidateExpSetup(ExpSetup,0); %from ScouseTom Repo

save([filepath 'Adult2016_TankCompSettings'],'ExpSetup');

%change to valid settings at top of BioSemi Range:

ExpSetup.Amp=240;
ExpSetup.Freq=1720;
%validate the new settings
[okflag,ExpSetup]=ScouseTom_ValidateExpSetup(ExpSetup,0); %from ScouseTom Repo

save([filepath 'AdultTank2016_ImagingSettings'],'ExpSetup');


%% Conductivities

disp('Setting Conductivities');
%conducivity for tank is saline everywhere, and then skull- there is no
%scalp or csf layers

Mesh.Nodes=1000*Mesh.Nodes;
ind_skull=find(Mesh.mat_ref==2);

%this
cnts=(Mesh.Nodes(Mesh.Tetra(:,1),:)+Mesh.Nodes(Mesh.Tetra(:,2),:)+Mesh.Nodes(Mesh.Tetra(:,3),:)+Mesh.Nodes(Mesh.Tetra(:,4),:))/4;
cnts=cnts(:,[1,3,2]);
origin=mean(cnts(ind_skull,:));

cnts=cnts-repmat(origin,length(cnts),1);

sigma=repmat(0.4,length(cnts),1);  %0.2% saline - closest from literature
sigma(Mesh.mat_ref==3)=0.00001; % make supports essentially infinite resistance
sigma0=0.0069*2; % reference conductivity

% adjust conductivity in X and Y 
sigskull=sigma0+abs(cnts(ind_skull,2))*0.87*sigma0/max(abs(cnts(ind_skull,2)))-abs(cnts(ind_skull,1))*0.32*sigma0/max(abs(cnts(ind_skull,1)));
% stick in full sigma vectors
sigma(ind_skull)=sigma0+abs(cnts(ind_skull,2))*0.87*sigma0/max(abs(cnts(ind_skull,2)))-abs(cnts(ind_skull,1))*0.32*sigma0/max(abs(cnts(ind_skull,1)));

disp('Writing VTKs');
% write VTK of fine mesh with matref and conductivities
writeVTKcell([filepath 'AdultFineMesh_MatRef'],Mesh.Tetra,Mesh.Nodes,Mesh.mat_ref);
writeVTKcell([filepath 'AdultFineMesh_Sigma'],Mesh.Tetra,Mesh.Nodes,sigma);

%% Write to dune format

disp('Writing Dune');
% dune exporter needs meters
dune_exporter(Mesh.Nodes(:,1:3)/1000,Mesh.Tetra(:,1:4),sigma,filepath,'Adult_Tank_2016.dgf',pos/1000,gnd_pos/1000,0);
