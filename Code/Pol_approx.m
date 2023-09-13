
% =========================================================================
% INPUT: 
% -----
% - Vg : strutcure that encode the function we want to approximate.
% 
%        * if flag = 0, Vg is a symbolic function. 
%        * if flag = 1, Vg contains Vg.basis (cell of handle functions) and 
%          Vg.coeff (Coefficient). For instance, 
%
%          f(x) = exp(x(1))-2*cos(x(2)) is encoded by 
%
%               o Vg.basis = {@(x) exp(x(1)),@(x) cos{x(2)}},
%               o Vg.coeff = [1,-2].
%          
% 
% - a : domain on which the polynomial approximation is computed. 
%       That is, the domain is [-a,a]^(dim).    
%
% - choice : Matrix of approximate eigenvectors of the Koopman operator
% 
% - order : order of the polynomial approximation of polynomial encoded by
%           Vg.
%
% - dim : dimension of the dynamical system
%
%
% OUTPUT: 
% ------
%
% - Expo contains the exponant of monomials of the polynomial approximation
%
% - Coeff contains the associated coefficients. 
% 
% - h is the bound on the error of approx: 
%
%       * Using Remez, it is estimated by h (the solution of the 
%         Lp program). 
%
%       * Using Taylor, it is not given but assumed to be given by the user
%         (see main.m) 
% 
% -------------------------------------------------------------------------

function [Coeff,Expo,h] = Pol_approx(Vg,a,choice,order,dim,flag)
    
    x = sym('x',[1 dim]);
    Expo = Exp_mon(dim,order);
    z = prod(x.^Expo,2);
    n = length(z);
    
    if(strcmp(choice,'minimax'))
        
        [alpha,h] = Remez_alg(Vg,n,z,a,flag,dim);   
        Coeff = value(alpha); 
    
%         indx2 = find(abs(Coeff)<1e-04); 
%         Coeff(indx2) = []; Expo(indx2,:) = []; 
%         indx2 = find(sum(Expo,2)==0); 
%         Coeff(indx2) = []; Expo(indx2,:) = [];

    else

        temp = taylor(Vg,x,zeros(dim,1),'Order',order+1);
        [coefficients, powers] = coeffs(temp,x);
        for i = 1:dim
            I = sym(ones(1,dim));
            I(i) = x(i);
            z = subs(powers,x,I);
            syms y positive
            Expo(:,i) = double(simplify(subs(log(z)./log(x(i)),x(i),y)));
        end
        Coeff = double(coefficients);
        h = []; 

    end

end