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
                    if (ang>360) % angle should be <45deg
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