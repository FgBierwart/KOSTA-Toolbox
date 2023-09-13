
% Assume E and C are monomials exponent and associated coefficient of the 
% polynomial ∑i v(i)*f(i). The following function compute its square. That
% is (∑i v(i)*f(i))^2 = ∑v(i)^2*f(i)^2 + 2∑_{i<j}v(i)*v(j)*f(i)*f(j)

function [Coeff,Expo] = pol_square_mon(C,E,dim)

    S1 = C.^2;
    T1 = 2*E;
    
    % computation of 2*Σ_{i<j}v_iv_j*f(i)*f(j)
    
    clear coeff_carr2 vec;
    s = 1;
    
    for i = 1:length(E)
        for j = (i+1):length(E)
            T2(s,:) = E(i,:)+E(j,:);
            S2(s) = 2*C(i)*C(j);
            s=s+1;
        end
    end
    
    B = unique(T2,'rows');
    D = zeros(length(B),1);
    
    for i = 1:size(B,1)
        h = find(sum(T2==B(i,:),2)==dim);
        D(i) = sum(S2(h));
    end
    
    % We are taking common elements to reduce the vector of coefficient 
    
    [temp,a,b] = intersect(T1,B,'rows');
    T1(a,:) = [];
    B(b,:) = [];
    vec = S1(a)+D(b);
    S1(a) = [];
    D(b) = [];
    
    % We regroup everything. Expo and Coeff contains new coefficients and
    % exponents of the new polynomial. 
    
    Expo = [T1;B;temp];
    Coeff = [S1;D;vec];

end