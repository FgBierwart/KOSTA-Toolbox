
function [Coeff,Expo,C_nrm_grad,E_nrm_grad] = Grad_fct_lyap(Coeff,Expo)

% According to the polynomial describe by Expo and Coeff, we are computing 
% the square of norm of the gradient of the polynomial. C_nrm_grad and 
% E_nrm_grad contain the final results 

dim = size(Expo,2); 

[Expo_g,Coeff_g] = Gradient_pol_mon(Expo,Coeff);
[C_nrm_grad,E_nrm_grad] = pol_square_mon(Coeff_g{1},Expo_g{1},dim);

    for j = 2:dim
        clear S T temp vec 
        [S,T] = pol_square_mon(Coeff_g{j},Expo_g{j},dim);
        [temp,a,x] = intersect(E_nrm_grad,T,'rows');
        E_nrm_grad(a,:) = [];
        T(x,:) = [];
        vec = C_nrm_grad(a)+S(x);
        C_nrm_grad(a) = [];
        S(x) = [];
        E_nrm_grad = [E_nrm_grad;T;temp];
        C_nrm_grad = [C_nrm_grad;S;vec]; 
    end

end