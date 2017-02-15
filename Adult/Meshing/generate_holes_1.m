function [electrodes, normals] = generate_holes_1 (Mesh)
%% Mesh has to have faces defined! use dubs2d for doing that

%this bit doesnt work kirill:

% if (~exist('Mesh.Faces') | ~exist('Mesh.Node_faces'))
%  ss=dubs3_2(Mesh.Tetra(:,1:4));
%  [Mesh.Node_faces,Mesh.Faces]=removeisolatednode_tri(Mesh.Nodes(:,1:3),ss);
%end
%% now it actually starts

vtx=Mesh.Nodes_faces; % NodeS faces - extra s
% vtx=vtx(:,[1,3,2]);

srf = Mesh.Faces;
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;
origin=mean(cnts);
%%

% THIS IS A CHEAT: we only need a half! but it will work for 4 quaters as
% well if you define everything correctly
srf=srf(cnts(:,1)>=0,:);
cnts=cnts(cnts(:,1)>=0,:);

% params
hole_distance=5;
hole_diam=1.1;

% eps=hole_diam is an option! this controls the epsilon-region
eps=mean((sum((vtx(srf(:,1),:)-vtx(srf(:,2),:)).^2,2)).^0.5)*0.8;

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
[el_s, n_side, sel] = make_stripe(first, 1, 2, sel, eps, hole_distance); % stripe of holes along x
electrodes=[electrodes;el_s];

% Past for the first stripe parallel to y
p=find(cnts(:,1)>=x0-eps & cnts(:,1)<=x0+eps);
sel=cnts(p,:);

% remove elements from queue
dist=(sum((sel-repmat(first,length(sel),1)).^2,2)).^0.5;
p=find(dist<4*eps);
sel(p,:)=[];

[el_f, n_front, sel] = make_stripe(first, 2, 1, sel, eps, hole_distance); % stripe of holes along y forward
electrodes=[electrodes;el_f];

[el_b, n_back, sel] = make_stripe(first, -2, 1, sel, eps, hole_distance); % stripe of holes along y backward
electrodes=[electrodes;el_b];

scatter3(electrodes(:,1),electrodes(:,2),electrodes(:,3),'b');daspect([1 1 1]);

%% now filling the quaters

% quater 1
sel=cnts(cnts(:,2)>y0,:);
[el,n] = fill_quater (el_f, el_s, sel, first, eps*2.2, hole_distance);
electrodes=[electrodes;el];

% quater 2
sel=cnts(cnts(:,2)<y0,:);
[el,n] = fill_quater (el_s, el_b, sel, first, eps*2, hole_distance);
electrodes=[electrodes;el];

% reflect the others as symmetry and cheese
el_flip=electrodes;
el_flip(:,1)=-el_flip(:,1);
el_flip(el_flip(:,1)>=-eps,:)=[];
electrodes=[electrodes;el_flip];


%% calculate normals
srf = Mesh.Faces;
cnts = (vtx(srf(:,1),:) + vtx(srf(:,2),:) + vtx(srf(:,3),:))/3;

%get two direction vectors X and Y 
t1=vtx(srf(:,1),:) - vtx(srf(:,2),:);
t2=vtx(srf(:,3),:) - vtx(srf(:,2),:);

%normal of surface given by cross product
norm_srf=cross(t1,t2);

%find the angle of the normal compared to the 
origins=cnts-repmat(origin, size(cnts,1),1);

%angle between normal and origin vector
cos_a=dot(norm_srf,origins,2)./(sum(norm_srf.^2,2).*sum(origins.^2,2)).^0.5;
angles=abs(acosd(cos_a));

%set all angles such as they all point outwards
norm_srf(angles>70,:)=-norm_srf(angles>70,:);

norm_srf=norm_srf./repmat(sum(norm_srf.^2,2),1,3);


%normals for each electrode/hole are mean of all faces close to the
%electrode/hole
kk=[]; 
for i=1:size(electrodes,1)
    %find surfaces close to hole
    dist=(sum((cnts-repmat(electrodes(i,:),length(cnts),1)).^2,2)).^0.5;
    p=find(dist<=4*eps);
    %sum the normals 
    normals(i,:)=sum(norm_srf(p,:),1);
    
    %reject any small or empty normals 
    if sum(normals(i,:).^2,2).^0.5 <=1e-10 || isempty (p)
        kk=[kk,i]; %add to remove idx
    end
end

electrodes(kk,:)=[];
normals(kk,:)=[];

%normalise the normals to a unit vector
normals=normals./repmat(sum(normals.^2,2).^0.5,1,3);

%chose a subset of these electrodes/holes - or dont
%p=find(electrodes(:,2)<-60);
p=1:size(electrodes,1);

%here we could check for weird curvature - if drastically different from
%nearest ones

figure;
scatter3(electrodes(p,1),electrodes(p,2),electrodes(p,3));
hold on;
quiver3 (electrodes(p,1),electrodes(p,2),electrodes(p,3),normals(p,1),normals(p,2),normals(p,3));
daspect([1 1 1]);


%% calculate hole size 

%this is the baseline diameter 1.1mm 
%and 4 holes per 10 x 10 mm square
%and 4 holes 
d=hole_diam; 


%the magic function to create the conductivities - from the impedance at
%the centre top of the skull - conductivity increases by 37% in x
%direction, and decreases by 15% in y direction (toward temporal)
diams=d+abs(electrodes(:,2))*0.37*d/max(abs(electrodes(:,2)))-abs(electrodes(:,1))*0.15*d/max(abs(electrodes(:,1)));

%% save SCAD file with cyllinders in mm
% writes a file to create all these cylinders for OpenSCAD - the boolean
% operation never really worked as the stls of the skulls had to be
% perfect, and they never were
%
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

%this creates the macro to run in FreeCAD to generate all the cylinders
%with the correct width, aligned with the normal of each point


%rearrange the ordering of the normals etc. to put it in the same order as
%freecad
O=[0,0,1];
normals=normals(:,[1,3,2]);
electrodes=electrodes(:,[1,3,2]);

%calculate the 
origins=repmat(O, size(normals,1),1);
cos_a=dot(normals,origins,2)./(sum(normals.^2,2).*sum(origins.^2,2)).^0.5;
axes=cross(origins,normals);
angles1=acosd(cos_a);
axes=axes./repmat(sum(axes.^2,2).^0.5,1,3);
radius=diams/2;
cyl_height='20.00';
F=fopen('Cyllinders_FreeCad_new.FCMacro','w');

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
        'App.ActiveDocument.C' num2str(i,'%4d') '.Radius=' num2str(radius(i),'%10.6f') '\n' ...
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
% radius=diams/2;
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
function [electrodes, n, sel] = make_stripe(first, dim_along, dim_cross, sel, eps, hole_distance)
% dimension- along which dimension, sign - left or right with right being
% positive;
% first-first electrode, not included in the output
% sel - selection of srf as a stripe
electrodes(1,:)=first;
n=2;
dd=abs(dim_along);

next=first;
status=0;
while (status~=2)
    dist=(sum((sel-repmat(next,size(sel,1),1)).^2,2)).^0.5;
    
    if (status==0) % place the first electrode in the direction
        if dim_along>0 % direction of placement
            p=find(dist>hole_distance-eps & dist<hole_distance+eps & sel(:,dd)>next(dd));
            p1=find(dist>hole_distance-3*eps & dist<hole_distance+3*eps & sel(:,dd)>next(dd));
        else
            p=find(dist>hole_distance-eps & dist<hole_distance+eps & sel(:,dd)<next(dd));
            p1=find(dist>hole_distance-3*eps & dist<hole_distance+3*eps & sel(:,dd)<next(dd));
        end
        status=1;
    else % all other electrodes
        p=find(dist>hole_distance-eps & dist<hole_distance+eps);
        p1=find(dist>hole_distance-3*eps & dist<hole_distance+3*eps);
    end
    if (~isnan(p)) % stop when could not find any elements left
        next=sel(p,:);
        sel(p1,:)=[];
        next=mean(next,1);
        next(dim_cross)=first(dim_cross);
        %check if the sharp edge
        place=1;
        if n>2
            v_p=electrodes(n-1,:)-electrodes(n-2,:);
            v_c=next-electrodes(n-1,:);
            cosa=dot(v_p,v_c)/(norm(v_p)*norm(v_c));
            ang=acosd(cosa);
            if (ang>30) % angle should be <45deg
                place=0;
            end
        end
        
        if place
            electrodes(n,:)=next;
            n=n+1;
        else
            status=2;
        end
    else
        status=2;
    end
end
electrodes(1,:)=[]; % we do not need first here
end

%%
function [electrodes, n] = fill_quater (el_d1, el_d2, sel, first, eps, hole_distance)
% this filles the quater between el_d1 and el_d2, with sel being the
% elements, first is our centroid
all_e=[first; el_d1; el_d2];
electrodes=[];
for i=1:size(all_e,1)
    temp=all_e(i,:);
    dist=(sum((sel-repmat(temp,length(sel),1)).^2,2)).^0.5;
    p=find(dist<2*eps);
    sel(p,:)=[];
end

el_raw=el_d2;

for i=1:size(el_d1,1)
    el_raw_temp=el_d1(i,:);
    for j=1:size(el_raw,1)
        
        dist1=(sum((sel-repmat(el_raw_temp(end,:),size(sel,1),1)).^2,2)).^0.5;
        dist2=(sum((sel-repmat(el_raw(j,:),size(sel,1),1)).^2,2)).^0.5;
        
        if (j==size(el_raw,1))
            p=find(dist1>hole_distance-eps & dist1<hole_distance+eps &dist2>=hole_distance-eps & dist2<=hole_distance+eps );
        else
            dist3=(sum((sel-repmat(el_raw(j+1,:),size(sel,1),1)).^2,2)).^0.5;
            p=find(dist1>hole_distance-eps & dist1<hole_distance+eps &dist2>=hole_distance-eps & dist2<=hole_distance+eps ...
                & dist3<=sqrt(2)*(hole_distance+eps) & dist3>=sqrt(2)*(hole_distance-eps)  );
        end
        
        
        if (~isnan(p))
            next=sel(p,:);
            next=mean(next,1);
            
            place=1;
            if size(el_raw_temp,1)>=2
                v_p=el_raw_temp(end,:)-el_raw_temp(end-1,:);
                v_c=next-el_raw_temp(end,:);
                cosa=dot(v_p,v_c)/(norm(v_p)*norm(v_c));
                ang=acosd(cosa);
                if (ang>45) % angle should be <45deg
                    place=0;
                end
            end
            
            if place
                el_raw_temp=[el_raw_temp;next];
                dist=(sum((sel-repmat(next,size(sel,1),1)).^2,2)).^0.5;
                sel(dist<=2*eps,:)=[];
            else
                j=size(el_raw,1);
            end
        else
            j=size(el_raw,1);
        end
        
    end
    electrodes=[electrodes;el_raw_temp(2:end,:)];
    el_raw=el_raw_temp(2:end,:);
end
n=size(electrodes,1);
end

%%
function Q = quat(vect,theta)

%Vector taken in [i j k], angle in degrees
%quaterions given in Format [w xi yj zk] or [w x y z]
%
%Example of using function:  ccc=quat([.3 .1 .6],45)
%Answer should be ->
%0.962730w   0.119633i   0.039878j   0.239266k

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
