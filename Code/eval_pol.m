
% This function evaluate the polynomial defined by Coeff and Expo over the
% points in v (v is a matrix where rows areassociated to the coordinates). 

function P = eval_pol(Coeff,Expo,v)

    P = zeros(size(v,1),1);
    
    for j = 1:length(P)
        P(j) = Coeff'*prod(v(j,:).^Expo,2);
    end

end