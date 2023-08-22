
% -------------------------------------------------------------------------
% This function computed the polynomial (or order 'order' approximation of 
% any handle function Vg over [-a,a]^2. Expo contains the exponant of 
% monomials and Coeff the associated coefficients. Here, h is the bound on 
% the error of approx: 
% 
%   * Using Remez, it is estimated by h (the solution of the Lp program). 
%   * Using Taylor, it is not given assumed to be given by the user
%     (see main.m) 
% 
% -------------------------------------------------------------------------

function [Coeff,Expo,h] = Pol_approx(Vg,a,choice,order,dim)

    syms x1 x2; x = [x1;x2]; 
    Expo = Exp_mon(dim,order);
    z = prod(transpose(x).^Expo,2);
    n = length(z);

    if(strcmp(choice,'minimax'))
        
        [alpha,h] = Remez_alg(Vg,n,z,a);   
        Coeff = value(alpha); 
    
%         indx2 = find(abs(Coeff)<1e-04); 
%         Coeff(indx2) = []; Expo(indx2,:) = []; 
%         indx2 = find(sum(Expo,2)==0); 
%         Coeff(indx2) = []; Expo(indx2,:) = [];

    else

        temp = taylor(Vg,x,zeros(dim,1),'Order',order+1);
        x = sym('x',[1 dim]);
        [coefficients, powers] = coeffs(temp,x);
        for i = 1:dim
            I = sym(ones(1,dim));
            I(i) = x(i);
            z = subs(powers,x,I);
            syms y positive
            Expo(:,i) = double(simplify(subs(log(z)./log(x(i)),x(i),y)));
        end
        Coeff = double(coefficients);
        h = 0; 

    end

end