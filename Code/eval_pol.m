
% This function evaluate the polynomial defined by Coeff and Expo over the
% points in v (v is a matrix where rows are associated to the coordinates). 

function P = eval_pol(Coeff,Expo,v)

dim = size(Expo,2); 
tmp = v(:,1)'.^Expo(:,1);

for j = 2:dim
    N = v(:,j)'.^Expo(:,j);
    tmp = tmp.*N; 
end

P = tmp'*Coeff;

end