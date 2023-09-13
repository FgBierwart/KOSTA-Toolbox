
% =========================================================================
% INPUT: 
% -----
% - basis : string caracter that contains the basis.
% 
% - s : structure that contains all information about the basis of
%       interest. E.g : 
%
%       * For monomials: s contains the order of monomials (Expo)
%       * For gaussian: s contains each basis functions encoded in s.f. 
%
% - Vec : Matrix of approximate eigenvectors of the Koopman operator
% 
% - indx : column of Vec that indicates the eigenvectors with closest 
%          eigenvalues of the one of jacobian matrix.
%
% - F_exp : Cell encoding the polynomial vector field or its approximation.
%           Each cell is a matrix associated to a componenent of the
%           vector. Each rows contains monomials and associated
%           coeffificient. E.g, F = [x(2) -2*x(1)+(1/3)*x(1)^3-x(2)];
% 
%           F_exp{1} = 0 1 1
%           F_exp{2} = |0 1 -1 |
%                      |1 0 -2 |
%                      |3 0 1/3|
%
%
% OUTPUT: 
% ------
% 
% - final_grid : Matrix where each rows is coordinates in X where the time
%                derivative of V is striclty negative. For each row, firsts 
%                column are the coordinates and the last one is the 
%                stepsize of the cell. E.g, [-0.5 -0.5 1]. Cell defined by 
%                the coordinates (-0.5,0.5), (-0.5,0.5), (0.5,-0.5),
%                (0.5,0.5). 
% 
% - grille_dot : Matrix with construction as final_grid contains the 
%                remaining cells in X which are not in final_grid.  
%
% - Expo : Matrix containing the exponents of monomials in V. 
% 
% - Coeff : Vector conntaining the coefficients of V associated to
%           monomials defined in Expo. 
%
% - Coeff_dotv : Same as Coeff for the time derivative of V. 
% 
% - Var_dotv : Same as Expo for the time derivative of V.
%
% =========================================================================

function [final_grid,grille_dot,Expo,Coeff,Coeff_dotv,Var_dotv] = AGM(Vec,indx,F_exp,basis,s,approx) 

    dim = length(indx);
    [Coeff_dotv,Var_dotv,C_grad_dot,E_grad_dot,Expo,Coeff,~,~] = Coeff_grad_dotv(Vec,indx,F_exp,basis,s,approx); 
    
    combs = Exp_mon(dim,dim-1); [a,~] = find(combs>1); combs(a,:)=[];
    combs = [combs;ones(1,dim)]; 
    
    % Creation of the grid.
    
    a = -1;
    Grille = Liste([a*ones(1,dim) 2],2); 
    
    nb_iter = 1;
    final_grid = []; grille_dot = []; 

    % Adaptive refinement of the grid.
    
    while(true)
    
    vec = []; div = 2; 
    
    A = repelem(Grille(:,1:end-1),2^dim,1);
    B = repmat(combs,size(Grille,1),1);
    C = repelem(Grille(:,end),2^dim,1); 
    
    % Here, tmp is a matrix so that any block of 2^n row corresponds to a 
    % cell and encodes (over its column) the value R_(r1,r2,...,rn) with 
    % (r1,...,rn) = (1,...,1), for each corner of the cell.  

    tmp = eval_pol(Coeff_dotv{1},Var_dotv{1},A+(B.*C));
    tmp = reshape(tmp,2^dim,[]).'; 
    
    % Here, we check cells that satisfies R_(r1,...,rn)(x) < 0 as well as
    % R_(r1,...,rn)(x) > 0, where x are corners of the cells.  

    ind1 = find(sum(sign(tmp),2)==(2^dim));
    ind2 = find(sum(sign(tmp),2)==-(2^dim)); 
    
    Cr = zeros(length(ind2),length(Coeff_dotv)); 
    
    if(~isempty(ind2))
        Cr(:,1) = min(-tmp(ind2,:),[],2); 
    end
    
    % Here, Cr contains at first the cells of the grid that satisfies one 
    % of the 2^n inequality.  

    vec = [vec;ind1]; 
    
    if(~isempty(ind2))
    
        l = 2; 
    
        while(l<=length(Coeff_dotv))
            
            tmp = eval_pol(Coeff_dotv{l},Var_dotv{l},A+(B.*C));
            tmp = reshape(tmp,2^dim,[]).'; 
    
            tmp2 = find(sum(sign(tmp),2)==(2^dim));
            vec = [vec;tmp2]; 
            
            t = find(sum(sign(tmp),2)==-(2^dim));
            
            [ind2,~,b] = intersect(t,ind2); 
            a = setdiff(1:size(Cr,1),b); 
            Cr(b,l) = min(-tmp(ind2,:),[],2);
            Cr(a,:) = [];  
            l = l+1; 
      
        end

        % At the end of the iteration, Cr contains all the cells satisfying 
        % the 2^n inequality. 
    
    end

    % For each of those cells, we check the criterion. 
    
    if(~isempty(ind2))    
         K = zeros(length(ind2),length(Coeff_dotv)); 
        for l = 1:size(K,2)
            K(:,l) = Bound_Grad_Pol(E_grad_dot{l},C_grad_dot{l},Grille(ind2,:),'max',combs);
            K(:,l) = (sqrt(K(:,l))*sqrt(dim)/2).*Grille(ind2,dim+1);
        end           
        q = find(sum(K<Cr,2)==size(K,2)); 
        grille_dot = [grille_dot;Grille(vec,:)]; 
        vec = [vec;ind2(q)];
        final_grid = [final_grid;Grille(ind2(q),1:dim+1)];
    end

    % We remove cells satisfying the 2^n criterion and cells where at least 
    % R_(r1,...,rn)(x) > 0 for one set of parameter (r1,...,rn) where x are 
    % corners of the cells (see ind1 and tmp2). This choice is motivated to 
    % avoid excessive division. 

    Grille(vec,:) = []; 
    
    % We stop the adaptive refinement where a fixed stepsize is reached. 

    if(~isempty(final_grid) && min(final_grid(:,end))<0.01) 
        grille_dot = [grille_dot;Grille]; 
        break;
    else
        Grille = Liste(Grille,div);
    end
    
    nb_iter=nb_iter+1;
    
    end

end
