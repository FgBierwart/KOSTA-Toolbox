
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
% - approx : structure that contains: approx.flag which is 0
%            for polynomial field or 1 for nonpolynomial vector fields, 
%            approx.err encodes the polynomial error (for non-polynomial 
%            vector field) in a same way as F_exp for each componenent of 
%            the vector field (E.g with taylor: err = c(i)||x||^d.) For 
%            polynomial vector fields, error is 0 for each component. It
%            also contains approx.basis which is 1 if the gradient of the
%            basis functions are approximated and 0 if there are given by
%            the user. 
%
%
% OUTPUT: 
% ------
%
% Asumme the general polynomial P where each monomials are recorder
% according the matrix 'E'. Each row contain exponant of monomials. The
% matrix 'C' are associated coefficient. Then: 
% 
% - Coeff_dotv : Coefficient of monomials in dV/dt computing with the 
%                approximate vector field F_exp 
% 
% - Var_dotv : Exponant of monomials in dV/dt computing with the approximate 
%              vector field F_exp 
% 
% - C_grad_dot : Coefficient of monomials in ||∇(dV/dT)||^2
% 
% - E_grad_dot : Exponent of monomials in ||∇(dV/dT)||^2
% 
% - Expo_dotg : Cell array where each row are the exponents of 
%               monomials of the components of ∇\dot{V}
% 
% - Coeff_dotg : Cell array where each components are the coefficient of 
%               monomials of the components of ∇\dot{V} 
% 
% - Expo_g : Cell array where each row are the exponents of monomials of 
%            the components of ∇V
% 
% - Coeff_g : Cell array where each components are the coefficient of 
%             monomials of the components of ∇V
% 
% - Expo : exponents of monomials in V (polynomial lyapunov candidate). 
%          When using gaussian functions, a minmax approximation a firstly 
%          computed. 
%
% - Coeff : Coefficient of monomials in V  
% 
% =========================================================================

function [Coeff_dotv,Var_dotv,C_grad_dot,E_grad_dot,Expo,Coeff,Expo_g,Coeff_g] = Coeff_grad_dotv(Vec,indx,F_exp,basis,s,approx)

dim = length(F_exp); 
choice = 'minimax'; a = 1; 

% Computation of the polynomial lyapunov function (exact with monomials or
% approximation using a minmax approximation for gaussian functions). 

if(~strcmp(basis,'monomials')) % minimax approx of the lyapunov function 
    
    order = 12;    
    Coeff = [real(Vec(:,indx)) imag(Vec(:,indx))]; 

    if(strcmp(choice,'minimax'))

        if(approx.basis~=0)

            % First possibility, we construct the Lyapunov approximation 
            % with handle function given in basis.  

            Vg = struct;
            Vg.basis = s.f;
            Vg.coeff = Coeff;
            flag = 1; 

        else

            % Or we use symbolic basis function to obtain the remez
            % version of the symbolic expression directly. (default) 
    
            g = s.fsym; 
            Vg = sum(sum(Coeff.*transpose(g),1).^2);
            flag = 0; 

        end

    elseif(strcmp(choice,'Taylor')) % if Taylor, we need symbolic. 
         g = s.fsym; 
         Vg.expr = sum(sum(Coeff.*transpose(g),1).^2); % symbolic
         flag = 0; 
    end
        
    [Coeff,Expo,~] = Pol_approx(Vg,a,choice,order,dim,flag);
    
    % tmp = find(abs(Coeff)<1e-03); 
    % Coeff(tmp) = []; Expo(tmp,:) = [];  

else
    [Coeff,Expo] = Lyapunov_pol(s.Expo,Vec,indx);
end

%  ============================= ETAPE 2 ================================

% Computation of ∇V(x) 

[Expo_g,Coeff_g] = Gradient_pol_mon(Expo,Coeff);

%  ============================ ETAPE 3 =================================

% Computation of ∇V(x).P(x) where P(x) is the polynomial approximation of 
% the vector field F(x).     

V = []; C = []; 

for k = 1:dim

    tmp = Expo_g{k};

    for i = 1:size(F_exp{k},1)    
        V = [V;tmp+F_exp{k}(i,1:end-1)];        
        C = [C;Coeff_g{k}.*F_exp{k}(i,end)];
    end

end

% if the vector field is not polynomial, we need to compute R_(r1,...,rn) 
% for some (r1,...,rn). the different r_i are encoded via combs.  

if(approx.flag==1) 

    N = 2^dim;
    combs = Exp_mon(dim,dim-1); [a,~] = find(combs>1); combs(a,:)=[];
    combs = [combs;ones(1,dim)]; combs(combs==0)=-1; combs = combs';

else 

    N = 1;

end
 

Var_dotv = cell(1,N); 
Coeff_dotv = cell(1,N); 

E_grad_dot = cell(1,N); 
C_grad_dot = cell(1,N); 

% if the vector field is not polynomial, Var_dotv{l} and Coeff_dotv{l} 
% encode the different R_(r1,...,rn) and C_grad_dot{l} and E_grad_dot{l}
% the associated gradient.

for l = 1:N
    
    Var_dotv{l} = V;
    Coeff_dotv{l} = C; 

    if(approx.flag==1) 

        for k = 1:dim
    
            tmp = Expo_g{k};
    
            for i = 1:size(approx.err.Expo{k},1)    
                Var_dotv{l} = [Var_dotv{l};tmp+approx.err.Expo{k}(i,:)];        
                Coeff_dotv{l} = [Coeff_dotv{l};Coeff_g{k}.*approx.err.Coeff{k}(i).*combs(k,l)];
            end
    
        end

    end

    H = unique(Var_dotv{l},'rows');
    H(:,dim+1) = zeros(size(H,1),1);
    
    for k = 1:size(H,1)
        id = ismember(Var_dotv{l},H(k,1:dim),'rows');
        H(k,end) = sum(Coeff_dotv{l}(id));   
    end
    
    Var_dotv{l} = H(:,1:dim);
    Coeff_dotv{l} = H(:,end);
    
    % Computation of ||∇R_(r1,...,rn)||^2. 
    
    [~,~,C_grad_dot{l},E_grad_dot{l}] = grad_norm(Var_dotv{l},Coeff_dotv{l}); 

end

end
