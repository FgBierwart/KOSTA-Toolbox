
% This function computes an upper bound of the polynomial (identified by 
% the matrices Expo and Coeff) over the cells in v. Here, v contains the 
% coordinates of the lower left vertices of the cells we are lookig at.  

function K = Bound_Grad_Pol(Expo,Coeff,v,type,combs)
    
    dim = size(Expo,2);
    vec = []; 

    m = size(v,1); 
    
    % Here, we check wich cells contains 0. 

    for i = 1:m
        for j = 1:dim
            tmp = find(v(i,j)<0 & v(i,j)+v(i,end)>0);
            if(~isempty(tmp))
                tmp = sort(find(Expo(:,j)>0)); 
                tmp = (tmp-1)*size(v,1)+i; 
                vec = [vec;tmp];
            end
        end
    end
    
    % Second: We construct the matrix of monomials for each cells. This 
    % contruction is similar to the one in AGM.m.         
        
    A = repelem(v(:,1:end-1),2^dim,1);
    B = repmat(combs,size(v,1),1);
    C = repelem(v(:,end),2^dim,1);
    
    v = A+(B.*C); 
    tmp = v(:,1)'.^Expo(:,1); 
        
    for j = 2:dim
        N = v(:,j)'.^Expo(:,j);
        tmp = tmp.*N;
    end
    
    tmp = reshape(tmp.',(2^dim),[]).'; 
    
    % For cells containing 0, the maximum (or minimum) over the cell is not
    % necessary at 0 so we have to add this quantity before computing the 
    % maximum (minimum).  
    
    tmp = [tmp tmp(:,1)];
    tmp(vec,end) = 0; 
    
    % Third: We compute the max or min for the bound. To do so, we identify 
    % cells where coefficients are negative and positive. For positive 
    % coefficient, we compute that maximum of the monomial and for negative 
    % coefficient we compute the minimum value of the monomial.      
    
    ta = find(Coeff>=0); i1 = []; 

    for i = 1:length(ta)
        i1 = [i1;(1:m)'+(ta(i)-1)*m]; 
    end
    tb = find(Coeff<0); i2 = []; 
    for i = 1:length(tb)
        i2 = [i2;(1:m)'+(tb(i)-1)*m];
    end

    t1 = max(tmp(i1,:),[],2); t2 = min(tmp(i2,:),[],2);
    
    t1 = reshape(t1,m,[]).';
    t2 = reshape(t2,m,[]).';
    
    K = t1'*Coeff(ta) + t2'*Coeff(tb);
    
    if(strcmp(type,'min'))
        K = -K;
    end    
     
end
            