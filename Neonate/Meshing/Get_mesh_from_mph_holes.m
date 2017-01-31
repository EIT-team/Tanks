files=dir('*.mph');
files={files.name};

for i=1:length(files)
    
model=mphload(files{i});
[stats,Data]=mphmeshstats(model);

Tetra=Data.elem{2};
Tetra=double(Tetra');
Mesh.Tetra=double(Tetra)+ones(size(Tetra,1),size(Tetra,2));
Mesh.Nodes=1000*Data.vertex';
Mesh.mat_ref=double(Data.elementity{2});
save([files{i}(1:end-3) '.mat'],'Mesh','-v7.3');

end