
% Isvalid takes a level set and ensure that it is a valid one. (see the
% paper for the explanation of how it is extracted).  

function [tmp,flag] = IsValid(G,V_grid,V,C,step,K,P)
    
    tmp = find(min(V_grid,[],2)<=V & max(V_grid,[],2)>=V);
    
    if(sum(G(tmp,end))~=length(tmp))
        
        flag = 0; 
        
    else
        
        vec = []; 
        
        % 1. Computation of the neighbours of the cells. 
        
        for i = 1:length(tmp)    
            vec = [vec; Neighbour(G,tmp(i))];
        end
        vec = unique(vec);          % remove duplicate cells 
        [~,b] = intersect(vec,tmp); % remove cells that cross the levelset. 
        vec(b) = []; 

        % 2. Regarding that we are looking for the top level set or bottom
        % level set, we only consider neighbour such either min > V or
        % max < V. 

        if(strcmp(step,'top'))
            vec = unique(vec(min(V_grid(vec,:),[],2)>V));
        else
            vec = unique(vec(max(V_grid(vec,:),[],2)<V));
        end

        % 2. Here we validate them according to the criterion.  

        flag = 1-CrossingCurve(G,vec,V_grid,V,C,step,K,P);


    end

end