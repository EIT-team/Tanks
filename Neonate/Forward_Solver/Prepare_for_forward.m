
%filepath=[filepath '/'];
filepath=[pwd filesep 'output' filesep];

disp('Loading Mesh');
%load mesh
load('Mesh_4p5_mln_forward.mat');
%load electrode positions
load('NNTank_Elec_pos.mat');

%% Ground node
disp('Finding Ground node');
%find ground node as node at bottom
gnd_pos=Mesh.Nodes(Mesh.Nodes(:,3)==min(Mesh.Nodes(:,3)),:);
gnd_pos=gnd_pos(1,1:3);

%% Protocol

disp('Finding Protocol');
%find "optimal" protocol
Prot_NN16=create_protocol(pos(1:end-1,:)); %from reconstruction repo

%find complete prt
[NN_Prt_full,NN_keep_idx,NN_rem_idx]=ScouseTom_data_findprt(Prot_NN16,32); %from ScouseTom Repo

save('NN2016Prt','NN_keep_idx','NN_Prt_full'); %for ref
save([filepath 'NN2016Prt'],'NN_keep_idx','NN_Prt_full');

%write protocol to textfile for forward solver
dlmwrite([filepath 'NN2016_Prt_full.txt'],NN_Prt_full);

%write protocol injections to textfile for ScouseTom System
dlmwrite([filepath 'NN2016_Prt_Injs.txt'],Prot_NN16);

%% ExpSetup

disp('Sorting ExpSetup');
%load previously used with with old protocol
load(['SkullTest2_log'],'ExpSetup');


%first make ExpSetup for Tank Comp Settings - this is at 300uA at 1khz. We
%cannot use this to make images, but we can use it to check tank is ok

%change protocol to new one
ExpSetup.Protocol=Prot_NN16;
%validate the new settings - ignore warning regarding too high current
[okflag,ExpSetup]=ScouseTom_ValidateExpSetup(ExpSetup,0); %from ScouseTom Repo

save([filepath 'NN2016_TankCompSettings'],'ExpSetup');

%change to valid settings at top of BioSemi Range:

ExpSetup.Amp=240;
ExpSetup.Freq=1720;
%validate the new settings
[okflag,ExpSetup]=ScouseTom_ValidateExpSetup(ExpSetup,0); %from ScouseTom Repo

save([filepath 'NN2016_ImagingSettings'],'ExpSetup');


%% Conductivities

disp('Setting Conductivities');
%conducivity for tank is saline everywhere, and then skull- there is no
%scalp or csf layers

sigma_back=0.4; %0.2% saline - closest from literature
sigma_sk=0.03; % from Pant2011

sigma_back_ref=1;
sigma_sk_ref=2;

sigma=repmat(sigma_back,length(Mesh.mat_ref),1); % saline
sigma(Mesh.mat_ref==sigma_sk_ref)= sigma_sk; % skull

disp('Writing VTKs');
% write VTK of fine mesh with matref and conductivities
writeVTKcell([filepath 'NNFineMesh_MatRef'],Mesh.Tetra,Mesh.Nodes,Mesh.mat_ref);
writeVTKcell([filepath 'NNFineMesh_Sigma'],Mesh.Tetra,Mesh.Nodes,sigma);



%% Write to dune format

disp('Writing Dune');
% dune exporter needs meters
dune_exporter(Mesh.Nodes(:,1:3)/1000,Mesh.Tetra(:,1:4),sigma,filepath,'NN2016_Tank.dgf',pos/1000,gnd_pos/1000,0);
