
% -------------------------------------------------------------------------
% This function convert a polynomial function into a matrix containing in
% firsts column the exponants of monomials and in the last column, the
% associated coefficients. 
% -------------------------------------------------------------------------

function F_exp = Convert_sym2mat(F_sym,dim)

F_exp = cell(1,length(F_sym)); 

for k = 1:length(F_sym)
    x = sym('x',[1 dim]);
    [coefficients, powers] = coeffs(F_sym(k),x);
    for i = 1:dim
        I = sym(ones(1,dim));
        I(i) = x(i);
        z = subs(powers,x,I);
        syms y positive
        F_exp{k}(:,i) = double(simplify(subs(log(z)./log(x(i)),x(i),y)));
    end
    F_exp{k}(:,dim+1) = double(coefficients);
end

end