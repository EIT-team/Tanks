
function [no,el]=removeisolatednode_tri(node,elem)

oid=1:size(node,1);       % old node index
 
ee=elem(:,1:3);
%ee=elem;
idx=setdiff(oid,ee(:)); % indices to the isolated nodes
idx=sort(idx);
delta=zeros(size(oid));   
delta(idx)=1;
delta=-cumsum(delta);     % calculate the new node index after removing the isolated nodes
oid=oid+delta;            % map to new index
el=oid(elem);             % element list in the new index
no=node;                  
no(idx,:)=[];             % remove the isolated nodes

end
