
% Given the the Grid (whose construction is similar to final_grid), this
% function computes all neighbor of the cells given by the row G(indx,:).

% This function is up to now implemented in 2D. 

function count = Neighbour(G,index)
    
    i = index;  
    
    count1 = find(G(:,1)+G(:,end-1)>=G(i,1) & G(:,1)<=G(i,1)+G(i,end-1)); 
    count2 = find(G(:,2)+G(:,end-1)>=G(i,2) & G(:,2)<=G(i,2)+G(i,end-1)); 
    
    a1 = find(G(count1,2)+G(count1,end-1)<=G(i,2) | G(count1,2)>=G(i,2)+G(i,end-1));
    
    count1(a1) = []; 
    
    a1 = find(G(count2,1)+G(count2,end-1)<=G(i,1) | G(count2,1)>=G(i,1)+G(i,end-1));
    count2(a1) = []; 
    
    count = unique([count1;count2]); 
    % count(count==i) = [];

end
