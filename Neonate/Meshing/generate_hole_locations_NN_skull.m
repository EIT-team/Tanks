function [electrodes, normals] = generate_hole_locations_NN_skull (Mesh)
%% Mesh has to have faces defined! use dubs2d for doing that

% if (~exist('Mesh.Faces') | ~exist('Mesh.Node_faces'))
%  ss=dubs3_2(Mesh.Tetra(:,1:4));
%  [Mesh.Node_faces,Mesh.Faces]=removeisolatednode_tri(Mesh.Nodes(:,1:3),ss);
% end
if (max(max(abs(Mesh.Node_faces)))<1)
    Mesh.Node_faces=1000*Mesh.Node_faces;
end
[Mesh.Node_faces,Mesh.Faces]=removeisolatednode_tri(Mesh.Node_faces,Mesh.Faces);
vtx=Mesh.Node_faces(:,1:3);
srf = Mesh.Faces;
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;
origin=mean(cnts);

vtx=vtx-repmat(origin,size(vtx,1),1);
cnts=cnts-repmat(origin,size(cnts,1),1);


eps=mean((sum((vtx(srf(:,1),:)-vtx(srf(:,2),:)).^2,2)).^0.5)*0.8;
srf=srf(cnts(:,3)>min(cnts(:,3))+eps,:);
cnts=cnts(cnts(:,3)>min(cnts(:,3))+eps,:);
%%
 TR=TriRep(srf, vtx(:,1),vtx(:,2),vtx(:,3));
 figure;
 h=trimesh(TR);
 set(h,'EdgeColor','k','FaceColor','none');
 daspect([1,1,1]);
 drawnow 
 hold on



% params
hole_distance=4;
hole_diam=1.5;

% eps=hole_diam is an option! this controls the epsilon-region



electrodes=[];

% here you can insert a code defining a centroid of the mesh! TODO
x0=0;
y0=0;
z0=max(cnts(:,3));

first=[x0,y0,z0]; % the first electrode/hole goes in the centroid

% Pass for first stripe parallel to x, we select a stripe parallel to x
% which containins the electrode
p=find(cnts(:,2)>=y0-eps & cnts(:,2)<=y0+eps);
sel=cnts(p,:);
dist=(sum((sel-repmat(first,length(sel),1)).^2,2)).^0.5;
[~,p]=min(dist);

first=sel(p,:); % the closest cnt to the first becomes first

% then we remove all cnts attached to the first from the queue so we do not
% consider them later
dist=(sum((sel-repmat(first,length(sel),1)).^2,2)).^0.5;
p=find(dist<4*eps);
sel(p,:)=[];

electrodes(1,:)=first; % yey, we've got first, now it is simple:
[el_s, n_side, sel] = make_stripe(first, 3, 2, sel, eps, hole_distance); % stripe of holes along x
electrodes=[electrodes;el_s];
    
[el_os, n_oside, sel] = make_stripe(first, -3, 2, sel, eps, hole_distance); % stripe of holes along x backward
electrodes=[electrodes;el_os];

    

 % Past for the first stripe parallel to y
p=find(cnts(:,1)>x0-eps & cnts(:,1)<=x0+eps);
sel=cnts(p,:);

% remove elements from queue
 dist=(sum((sel-repmat(first,length(sel),1)).^2,2)).^0.5;
 p=find(dist<3*eps);
 sel(p,:)=[];

[el_f, n_front, sel] = make_stripe(first, 2, 1, sel, eps, hole_distance); % stripe of holes along y forward
electrodes=[electrodes;el_f];

[el_b, n_back, sel] = make_stripe(first, -2, 1, sel, eps, hole_distance); % stripe of holes along y backward
electrodes=[electrodes;el_b];

scatter3(electrodes(:,1),electrodes(:,2),electrodes(:,3),'r','filled');
drawnow



%% now filling the quaters

% quater 1
sel=cnts(cnts(:,2)>y0,:);
[el,n] = fill_quater (el_s, el_f, sel, first, eps*3, hole_distance); 
electrodes=[electrodes;el];
scatter3(el(:,1),el(:,2),el(:,3),'b','filled');
 drawnow

% quater 2
sel=cnts(cnts(:,2)<y0,:);
[el,n] = fill_quater (el_s, el_b, sel, first, eps*3, hole_distance); 
electrodes=[electrodes;el];
scatter3(el(:,1),el(:,2),el(:,3),'b','filled');
 drawnow
sel=cnts(cnts(:,2)>y0,:);
[el,n] = fill_quater (el_os, el_f, sel, first, eps*3, hole_distance); 
electrodes=[electrodes;el];
scatter3(el(:,1),el(:,2),el(:,3),'b','filled');
 drawnow

sel=cnts(cnts(:,2)<y0,:);
[el,n] = fill_quater (el_os, el_b, sel, first, eps*3, hole_distance); 
electrodes=[electrodes;el];
scatter3(el(:,1),el(:,2),el(:,3),'b','filled');
 drawnow


%% calculate normals
vtx=vtx+repmat(origin,size(vtx,1),1);
electrodes=electrodes+repmat(origin,size(electrodes,1),1);

srf = Mesh.Faces;
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;

t1=vtx(srf(:,1),:) - vtx(srf(:,2),:);
t2=vtx(srf(:,3),:) - vtx(srf(:,2),:);

norm_srf=cross(t1,t2);

origins=cnts-repmat(origin, size(cnts,1),1);

cos_a=dot(norm_srf,origins,2)./(sum(norm_srf.^2,2).*sum(origins.^2,2)).^0.5;
angles=abs(acosd(cos_a));

norm_srf(angles>70,:)=-norm_srf(angles>70,:);

norm_srf=norm_srf./repmat(sum(norm_srf.^2,2),1,3);

kk=[];
for i=1:size(electrodes,1)
    dist=(sum((cnts-repmat(electrodes(i,:),length(cnts),1)).^2,2)).^0.6;
    p=find(dist<=4*eps);
    normals(i,:)=sum(norm_srf(p,:),1);
    if sum(normals(i,:).^2,2).^0.5 <=1e-10 || isempty (p)
        kk=[kk,i];
    end
end

electrodes(kk,:)=[];
normals(kk,:)=[];

normals=normals./repmat(sum(normals.^2,2).^0.5,1,3);

%p=find(electrodes(:,2)<-60);
p=1:size(electrodes,1);
figure;
scatter3(electrodes(p,1),electrodes(p,2),electrodes(p,3),'b','filled');
%hold on;
%quiver3 (electrodes(p,1),electrodes(p,2),electrodes(p,3),normals(p,1),normals(p,2),normals(p,3));
daspect([1 1 1]);

%% save SCAD file with cyllinders in mm
% O=[0,0,1];
% normals=normals(:,[1,3,2]);
% electrodes=electrodes(:,[1,3,2]);
% 
% origins=repmat(O, size(normals,1),1);
% cos_a=dot(normals,origins,2)./(sum(normals.^2,2).*sum(origins.^2,2)).^0.5;
% axes=cross(origins,normals);
% angles=acosd(cos_a);
% radius=hole_diam/2;
% cyl_height='10';
% F=fopen('Cyllinders_SCAD.scad','w');
% for i=1:size(electrodes,1)
%     translate_vector = ['[' num2str(electrodes(i,1),'%10.6f') ',' num2str(electrodes(i,2),'%10.6f') ',' num2str(electrodes(i,3),'%10.6f') ']'];
%     rotate_vector =  ['[' num2str(axes(i,1),'%10.6f') ',' num2str(axes(i,2),'%10.6f') ',' num2str(axes(i,3),'%10.6f') ']'];
%     
%     fprintf(F,['translate(' translate_vector ' ) {\n' ...
%              '   rotate(' num2str(angles(i),'%10.6f') ',' rotate_vector ') {\n' ...
%              '       cylinder(h=' cyl_height ',r=' num2str(radius,'%10.6f') ',$fn=25, center=true);\n' ...
%              '}}\n']);
% end
% fclose(F);

%% save FCMacro file with cyllinders in mm

disp('Writing FreeCAD Macro');
O=[0,0,1];
%normals=normals(:,[1,3,2]);
%electrodes=electrodes(:,[1,3,2]);

origins=repmat(O, size(normals,1),1);
cos_a=dot(normals,origins,2)./(sum(normals.^2,2).*sum(origins.^2,2)).^0.5;
axes=cross(origins,normals);
angles1=acosd(cos_a);
axes=axes./repmat(sum(axes.^2,2).^0.5,1,3);
radius=hole_diam/2;
cyl_height='15.00';
F=fopen('output/Cyllinders_FreeCad_new.FCMacro','w');

fprintf(F,['import FreeCAD\n' ...
          'import Part\n' ...
          'import Part,PartGui\n' ...
          'from FreeCAD import Base\n'] );

for i=1:size(electrodes,1)
    base=electrodes(i,:)-normals(i,:).*10;
    translate_vector = ['(' num2str(base(1),'%10.6f') ',' num2str(base(2),'%10.6f') ',' num2str(base(3),'%10.6f') ')'];
    Q1=quat(axes(i,:),angles1(i));
    rotate_vector1 =  ['(' num2str(Q1(2),'%10.6f') ',' num2str(Q1(3),'%10.6f') ',' num2str(Q1(4),'%10.6f') ',' num2str(Q1(1),'%10.6f') ')' ];
        
    fprintf(F,['App.ActiveDocument.addObject("Part::Cylinder","C' num2str(i,'%4d') '")\n' ...
               'App.ActiveDocument.C' num2str(i,'%4d') '.Radius=' num2str(radius,'%10.6f') '\n' ...
               'App.ActiveDocument.C' num2str(i,'%4d') '.Height=' cyl_height '\n' ...
               'App.ActiveDocument.C' num2str(i,'%4d') '.Angle=360.00 \n' ...
               'App.ActiveDocument.C' num2str(i,'%4d') '.Placement=Base.Placement(Base.Vector' translate_vector ',Base.Rotation' rotate_vector1 ')\n'] );
    
    
end

fprintf(F,'App.ActiveDocument.recompute()\n');

fclose(F);


%% save text file with cyllinders position-axis-angle  in mm
% O=[0,0,1];
% normals=normals(:,[1,3,2]);
% electrodes=electrodes(:,[1,3,2]);
% 
% origins=repmat(O, size(normals,1),1);
% cos_a=dot(normals,origins,2)./(sum(normals.^2,2).*sum(origins.^2,2)).^0.5;
% axes=cross(origins,normals);
% angles=acosd(cos_a);
% axes=axes./repmat(sum(axes.^2,2).^0.5,1,3);
% radius=hole_diam/2;
% cyl_height='10';
% F=fopen('Cyllinders.txt','w');
% for i=1:size(electrodes,1)
%     base=electrodes(i,:)-normals(i,:).*5;
%     translate_vector = [ num2str(base(1),'%10.6f') ',' num2str(base(2),'%10.6f') ',' num2str(base(3),'%10.6f') ];
%     rotate_vector =  [ num2str(axes(i,1),'%10.6f') ',' num2str(axes(i,2),'%10.6f') ',' num2str(axes(i,3),'%10.6f') ];
%     
%     fprintf(F,[translate_vector ',' rotate_vector ',' num2str(angles(i),'%10.6f') '\n']);
% end
% fclose(F);

end

%%

%% 


%%
%Vector taken in [i j k], angle in degrees
%quaterions given in Format [w xi yj zk] or [w x y z] 
%
%Example of using function:  ccc=quat([.3 .1 .6],45)  
%Answer should be ->
%0.962730w   0.119633i   0.039878j   0.239266k

function Q = quat(vect,theta)
	Q = zeros(1,4);
	xyz_and_angle=[vect theta] ;   %displays inputed data and theta angle before deg to rad conversion
	%thetaconvert = (theta*pi/180)./norm(theta*pi/180) %convert deg to rad
	thetaconvert=theta*pi/180;
	px=vect(1,1);py=vect(1,2);pz=vect(1,3); 
	Q(1,1) = cos(thetaconvert/2); %w value
	Q(1,2)= px*sin(thetaconvert/2); % x value
	Q(1,3)= py*sin(thetaconvert/2); % y value
	Q(1,4)= pz*sin(thetaconvert/2); % z value
	
	q1w=Q(1,1);%qw value
	q1x=Q(1,2);%qx value
	q1y=Q(1,3);%qy value
	q1z=Q(1,4);%qz value
	
	%Normalize Quaternion due to floating point errors
     	qnorm=sqrt(q1w^2+q1x^2+q1y^2+q1z^2) ;
     	Q=Q./qnorm ;%normalize Q array
end
