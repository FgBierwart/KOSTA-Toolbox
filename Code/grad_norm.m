
function [Coeff,Expo,C_grad_dot,E_grad_dot] = grad_norm(Expo,Coeff)

% According to the polynomial describe by Expo and Coeff, we are computing 
% the square of norm of the gradient of the polynomial. C_nrm_grad and 
% E_nrm_grad contain the final results 

dim = size(Expo,2); 
[Expo_dotg,Coeff_dotg] = Gradient_pol_mon(Expo,Coeff);
[C_grad_dot,E_grad_dot] = pol_square_mon(Coeff_dotg{1},Expo_dotg{1},dim);
       
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