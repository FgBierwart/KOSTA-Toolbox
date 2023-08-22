
% =========================================================================
% INPUT: 
% -----
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
% When using gaussian functions, a minmax approximation a firstly computed
%
% - Coeff : Coefficient of monomials in V  
% 
% =========================================================================

function [Coeff_dotv,Var_dotv,C_grad_dot,E_grad_dot,Expo,Coeff,Expo_dotg,Coeff_dotg,Expo_g,Coeff_g] = Coeff_grad_dotv(Vec,indx,F_exp,basis,s)

dim = length(indx); 
choice = 'minimax'; a = 1; 

% Computation of the polynomial lyapunov function (exact with monomials or
% approximation using a minmax approximation for gaussian functions). 

if(strcmp(basis,'gaussian')) 
    gamma = s.a; don = s.d'; order = 12;
    x = sym('x',[1 dim]);
    Coeff = [real(Vec(:,indx)) imag(Vec(:,indx))]; 
    g = exp(-gamma^2*sum((x-don).^2,2)); 
    Vg = sum(sum(Coeff.*g,1).^2); 
    [Coeff,Expo,~] = Pol_approx(Vg,a,choice,order,dim); 
else
    combs = Exp_mon(dim,s.a); combs(1,:) = [];
    [Coeff,Expo] = Lyapunov_pol(combs,Vec,indx);
end

% tmp = find(abs(Coeff)<1e-06); 
% Coeff(tmp) = []; Expo(tmp,:) = [];  

%  ============================= ETAPE 2 ================================

% Computation of ∇V(x) 

[Expo_g,Coeff_g] = Gradient_pol_mon(Expo,Coeff);

%  ============================ ETAPE 3 =================================

% CALCUL DE dV/dt 

% Var_dotv   : Contains monomials of Contient les exposants de dV/dt 
% Coeff_dotv : Contient les coefficients de dV/dt

% -------------------------------------------------------------------------

% We first compute F(i)*∇V(x)_i

Coeff_dotv = [];
Var_dotv = [];

for k = 1:dim
    tmp = Expo_g{k};
    for i = 1:size(F_exp{k},1)
        Var_dotv = [Var_dotv;tmp+F_exp{k}(i,1:end-1)];
        Coeff_dotv = [Coeff_dotv;Coeff_g{k}.*F_exp{k}(i,end)];
    end
end

% Computation of ||∇\dot{V}||^2

H = unique(Var_dotv,'rows');
H(:,dim+1) = zeros(size(H,1),1);

for k = 1:size(H,1)
    id = ismember(Var_dotv,H(k,1:dim),'rows');
    H(k,end) = sum(Coeff_dotv(id));   
end

Var_dotv = H(:,1:dim);
Coeff_dotv = H(:,end);

[Expo_dotg,Coeff_dotg] = Gradient_pol_mon(H(:,1:end-1),H(:,end));
[C_grad_dot,E_grad_dot] = pol_square_mon(Coeff_dotg{1},Expo_dotg{1},dim);

% Boucle servant à additionner les polynomes au carre (different par raport
% a V(x) car ici les variables ne sont pas les memes pour chaque
% composantes (∇V.F pour rappel) 

for j = 2:dim
    clear S T temp vec 
    [S,T] = pol_square_mon(Coeff_dotg{j},Expo_dotg{j},dim);
    [temp,a,x] = intersect(E_grad_dot,T,'rows');
    E_grad_dot(a,:) = [];
    T(x,:) = [];
    vec = C_grad_dot(a)+S(x);
    C_grad_dot(a) = [];
    S(x) = [];
    E_grad_dot = [E_grad_dot;T;temp];
    C_grad_dot = [C_grad_dot;S;vec]; 
end

end
