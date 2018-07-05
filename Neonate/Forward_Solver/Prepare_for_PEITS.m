% this is the high resolution mesh ~6m elements suitable for PEITS

savepath=[pwd filesep 'output' filesep];

disp('Loading Mesh');
Mesh=loadmesh('../Meshing/output/high_res/NNmesh_fine'); %load mesh

%in case there is extra electrode added by the mesher
Mesh.elec_pos=Mesh.elec_pos(1:33,:);


%% Protocol

disp('Finding Protocol');
%find "optimal" protocol
Prot_NN16=create_protocol(Mesh.elec_pos(1:end-1,:)); %from reconstruction repo

%find complete prt
[NN_Prt_full,NN_keep_idx,NN_rem_idx]=ScouseTom_data_findprt(Prot_NN16,32); %from ScouseTom Repo

% save('NN2016Prt','NN_keep_idx','NN_Prt_full'); %for ref
% save([savepath 'NN2016Prt'],'NN_keep_idx','NN_Prt_full');

%write protocol to textfile for forward solver
% dlmwrite([savepath 'NN2016_Prt_full.txt'],NN_Prt_full);

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

save([savepath 'NNPEITSMesh'],'Mesh');
