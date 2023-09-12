
function [Coeff_V,Expo_V] = Lyapunov_pol(combs,Vec,indx)

% Computation of the Lyapunov candidate V(x) = [Σ(phi_i^2)]^2 with 
% phi_i = Σ(v_i*f(i)) where f(i) are basis functions (given by combs) with 
% v_i the coefficient associated in Vec (see Eigenfunction.m) for more
% details. 

% 1. Computation of Σv_i²f(i)². 

dim = size(combs,2);
Real = real(Vec(:,indx));
Imag = imag(Vec(:,indx));

coeff_carr1 = sum(Real.^2,2)+sum(Imag.^2,2);
var_coeff1 = 2*combs;

% 2. Computation of calcul de 2*Σ_{i<j}v_iv_j*f(i)*f(j)

clear coeff_carr2 var_coeff2 vec;
s = 1;
n = length(combs);

for i = 1:n
    for j = (i+1):n
        var_coeff2(s,:) = combs(i,:)+combs(j,:);
        coeff_carr2(s) = 2*[Real(i,:) Imag(i,:)]*[Real(j,:) Imag(j,:)]';
        s=s+1;
    end
end

coeff_carr2 = coeff_carr2';

B = unique(var_coeff2,'rows');
D = zeros(length(B),1);

for i = 1:size(B,1)
    h = find(sum(var_coeff2==B(i,:),2)==dim);
    D(i) = sum(coeff_carr2(h));
end

% We are taking common elements to reduce the vector of coefficient 

[temp,a,b] = intersect(var_coeff1,B,'rows');
var_coeff1(a,:) = [];
B(b,:) = [];
vec = coeff_carr1(a)+D(b);
coeff_carr1(a) = [];
D(b) = [];

% We regroup everything. Expo_V and Coeff_V contain coefficient of V.

Expo_V = [var_coeff1;B;temp];
Coeff_V = [coeff_carr1;D;vec];

end