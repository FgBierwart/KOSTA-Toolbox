
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
% - indx : cell (index of G) for which we want to determine wheter 
%          or not it crosses the level curve V_star. Compare to crrosing
%          curve, indx here only contains one element which corresponds to
%          only one cell. 
%
% - C : array that contains Expo and Coeff, which are the exponent and 
%       coefficient of the candidate polynomial lyapunov function. See for 
%       documentation of Coeff_grad_dotv.m for more details.    
%
% - type : string character which indicates if the cells we are looking at 
%          need to be checked from above (min(V) > gam_2) or from below 
%          (max(V) < gam_1). Indeed, during the bisection, for a set of 
%          parameters gam_1 and gam_2, we check the cells that intersect 
%          the level set and the possible neigbours that might intersect 
%          the cell. For those neighbours, depending the step of bisection 
%          (bisection from 'top' or 'bottom', we only need to chek beighbor 
%          such that (max(V) < gam_1) or (min(V) > gam_2). 
%          
%          See IsValid.m for more details. 
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
% - flag : 0 if the level curve does not cross the cell and 1 otherwise.  
% 
% =========================================================================

function flag = CrossingCurve_pct(G,V_grid,indx,V_star,C,type,K,P)

Expo = C{1}; Coeff = C{2}; 

% Generating equidistinct points over hypercube. See the paper for details. 
N = 100; t = Hypercube_data(size(Expo,2),N,G(indx,:)',G(indx,end-1));
K = sqrt(K)/(2*(N-1))*G(indx,end-1); 

P = P(t); 

flag = 1; 

dim = 2; combs = Exp_mon(dim,dim-1); [a,~] = find(combs>1); combs(a,:)=[];
combs = [combs;ones(1,dim)];

if(strcmp(type,'top')) 

    % If we are in the case for a cell where min(V) > V_star, we compute
    % another bound (Kmin) on the minimum value of V over the cell (since 
    % it is a polynomial). Here, min(P)-K corresponds to a worst decreases 
    % scenario over the cell. Thus, if either min(P)-K > V_star or 
    % Kmin > V_star, we can ensure that the level curve does not cross the 
    % cell. 

    Kmin = Bound_Grad_Pol(Expo,-Coeff,G(indx,1:end-1),'min',combs);
    vec = min(P)-K; 

    if(min(V_grid(indx,:))>V_star && (vec>V_star || Kmin>V_star))
        flag = 0;
    end

elseif(strcmp(type,'bot'))

    % If we are in the case for a cell where max(V) < V_star, we compute
    % another bound (Kmax) on the maximum value of V over the cell (since 
    % it is a polynomial). Here, max(P)+K corresponds to a worst increases 
    % scenario over the cell. Thus, if either max(P)+K < V_star or 
    % Kmax < V_star, we can ensure that the level curve does not cross the 
    % cell.  
    
    vec = max(P)+K; 
    Kmax = Bound_Grad_Pol(Expo,Coeff,G(indx,1:end-1),'max',combs);

    if(max(V_grid(indx,:))<V_star && (vec<V_star || Kmax<V_star))
        flag = 0;  
    end

end

end
