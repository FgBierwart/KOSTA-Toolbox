
% =========================================================================
% INPUT: 
% -----
% - basis : string caracter that contains the basis.
% 
% - s : structure that contains all information about the basis of
%       interest. E.g : 
%       * For monomials : s contains the order of monomials
%       * For gaussian, s contains the value of gamma, the domain of the 
%         centers and the centers themselves.
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

function [final_grid,grille_dot,Expo,Coeff,Coeff_dotv,Var_dotv] = AGM(Vec,indx,F_exp,basis,s) 

% -- ETAPE 1 --
%
% Creation d'un vecteur d'entree 0/1 permettant de determiner les coins de
% la grille. 

dim = length(indx);
[Coeff_dotv,Var_dotv,C_grad_dot,E_grad_dot,Expo,Coeff] = Coeff_grad_dotv(Vec,indx,F_exp,basis,s);

combs = Exp_mon(dim,dim-1); [a,~] = find(combs>1); combs(a,:)=[];
combs = [combs;ones(1,dim)]; 

% -- ETAPE 2 -- 
% 
% 1. Creation de la grille initiale 

a = -1;
Grille = Liste([a*ones(1,dim) 2],2); 

nb_iter = 1;
final_grid = []; grille_dot = []; 

while(true)

    div = 2; % On divise chaque cellule par 2 
    vec = [];    

    % For each cells, we are going to compute the bound of the polynomial
    
    z = zeros(size(Grille,1),div^dim);

    for i = 1:size(Grille,1)
        v = Grille(i,1:end-1)+combs*Grille(i,end);
        z(i,:) = eval_pol(Coeff_dotv,Var_dotv,v); 
    end

    % Here, we are going to remove cells where we know that dot_v > 0. 
    
    ind = find(sum(sign(z),2)==2^dim); 
    vec = [vec;ind];

    % Here, we remove cells where dot_V is negative everywhere inside 

    %[row,~] = find(z==0);
    %row = sort(row);
    ind = find(sum(sign(z),2)==-2^dim);

    if(~isempty(ind))
        K = zeros(length(ind),1);
        for i = 1:length(ind)
            v = Grille(ind(i),1:dim)+combs*Grille(ind(i),dim+1);
            K(i) = Bound_Grad_Pol(E_grad_dot,C_grad_dot,v,Grille(ind(i),dim+1),'max');
            K(i) = sqrt(K(i))*sqrt(dim)/2*Grille(ind(i),dim+1); 
            % K(i) = Bound_Grad_Pol(Var_dotv,Coeff_dotv,v,Grille(ind(i),dim+1));
        end
        z = min(abs(z(ind,:)),[],2);
        q = find(K<z);
        % q = find(K<0); 
        grille_dot = [grille_dot;Grille(vec,:)]; 
        vec = [vec;ind(q)];

        final_grid = [final_grid;Grille(ind(q),1:dim+1)];
    end
    
    Grille(vec,:) = []; 
    if(~isempty(final_grid) && min(final_grid(:,end))<0.01)
        grille_dot = [grille_dot;Grille]; 
        break;
    else
        Grille = Liste(Grille,div);
    end
    %disp(nb_iter); 
    nb_iter=nb_iter+1;

end

end
