
% =========================================================================
% INPUT: 
% -----
% 
% - G : grid where the first block contain final_grid (see AGM.m) and the 
%       second one is grille_dot. The last coloumn of G is 1 for row 
%       associated to final_grid and 0 otherwise.    
%
% - V_grid : Matrix so that each row corresponds to a cell in G and columns 
%            are the values of the candidate lyapunov function V at each 
%            cells corner. For instance, if the cell of corrdinates are 
%            (0.5,0.5) with stepsize 0.1 is encoded by the 45th row of of 
%            G, then the 45th row of V is the value of V over the corners 
%            of this cell. 
%
% - C : see CrossingCurve_pct.m.   
%
% - step : string character which indicates if the cell we are looking at 
%          needs to be checked from above (min(V) > gam_2) or from below 
%          (max(V) < gam_1). Indeed, during the bisection, for a set of 
%          parameters gam_1 and gam_2, we check the cells that intersect 
%          the level set and the possible neigbours that might intersect 
%          the cell. For those neighbours, depending the step of bisection 
%          (bisection from 'top' or 'bottom', we only need to chek beighbor 
%          such that (max(V) < gam_1) or (min(V) > gam_2). Here, gam_2 =
%          gam_1 = V_star. We check if the leve curve is valid. 
%
%
% - K : This vector encodes bounds on the norm of the gradient of V over
%       the different cells in G. This is needed for the bisection.  
% 
% - P : handle function of the lyapunov candidate V
% 
% 
% OUTPUT: 
% ------
%
% - tmp : vector of index (cells) that cross the cell. More precisely,
%         these are cells suct that min(V) < V_star and max(V) > V_star. 
% 
% - flag : 1 if the level curve V_star only cross cells that are 'valid', 
%          that is cell in final_grid.   
% 
% =========================================================================

function [tmp,flag] = IsValid(G,V_grid,V_star,C,step,K,P)
    
    % Here, we check cells for which we are sure they cross the level
    % curve. These cells are such that min(V) < gam_1 and max(V) > gam_2.

    % Illustration of 4 cells with 1 cells that cross the level curve 
    %
    % *----*----* 
    % | 1  | 2  |       
    % *----*----*
    % | 3  | 4 _|_ 
    % *----*-/--*
    %       /
    %       | <--- V_star
    %
    % In this situation, V_star cross cell 4 since max(V) > V_star and
    % min(V) < V_star.

    tmp = find(min(V_grid,[],2)<=V_star & max(V_grid,[],2)>=V_star);
    
    if(sum(G(tmp,end))~=length(tmp)) 

        % The level set is not valid if at least one of the cell that cross 
        % the level set is not in final_grid.  
        
        flag = 0; 
        
    else
        
        % according to the criteria, max(V) > V_star and min(V) < V_star, 
        % we can miss cell that cross the cell. for instance, if we have : 

        % *------*------*
        % | 3 _  | 4____|_ 
        % *--|-|-*-/----*   
        %   /  |_/ 
        %   |               
        %
        % In this situation, 3 does not satisfy the criterion but 
        % crosses the level curve. Thus, in order to make sure that the 
        % level curve is valid, we will idnetify all cells such as 4 and 
        % all its neigbours that might cross the cells such as 3.   

        vec = []; 
        
        % Computation of the neighbours of the cells. 
        
        for i = 1:length(tmp)    
            vec = [vec; Neighbour(G,tmp(i))];
        end
        vec = unique(vec);          % remove duplicate cells 
        [~,b] = intersect(vec,tmp); % remove cells that cross the levelset. 
        vec(b) = []; 
        
        % Regarding that we are looking for the top level set (gam_2) or 
        % bottom level set (gam_1), we only consider neighbours such that 
        % either (min > V) or (max < V). Indeed, if try to identify gam_2, 
        % we only need to verify neigbours such that min(V) > V_star. For 
        % gam_1, that will be neighbours for which max(V) < V_star. The 
        % validation of the cells is done according to to Crossin_Curve.m. 
        %
        % Remark 1. We can only check among those particular neigbhours 
        % (regarding type = 'bot' or 'top') since during the bisection,
        % once two level curves gam_1 and gam_2 are identified and valid, 
        % we check all cells between these curves so that we are sure that 
        % the region defined by gam_1 < V < gam_2 is valid. In other 
        % worlds, we do not miss cell.
        %
        % Remark 2. For the reference level set, we need to verify all
        % neighbours since it will be the references. 

        if(strcmp(step,'top'))
            vec = unique(vec(min(V_grid(vec,:),[],2)>V_star));
        else
            vec = unique(vec(max(V_grid(vec,:),[],2)<V_star));
        end

        % Here we validate them according to the criterion.  

        flag = 1-CrossingCurve(G,vec,V_grid,V_star,C,step,K,P);

    end

end