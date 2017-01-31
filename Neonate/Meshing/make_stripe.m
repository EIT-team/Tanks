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
            if (ang>360) % angle should be <45deg
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
