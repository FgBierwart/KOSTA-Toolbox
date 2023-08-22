
% This function computes an upper of the polynomial (identified by the
% matrices Expo and Coeff) over the cell v. Here, v contains the differents 
% vertices of the cell we are lookig at. 

function K = Bound_Grad_Pol(Expo,Coeff,v,step,type)
    
    dim = size(Expo,2);

    % Here, we compute the values of the monomials at each corners. 
    tmp = v(:,1)'.^Expo(:,1); 
        
    for j = 2:dim
        N = v(:,j)'.^Expo(:,j);
        tmp = tmp.*N;
    end

    % Here, we check if the cell contains 0 for one variables. 
    l = zeros(dim,1);

    for i = 1:dim
        if(v(1,i)<0 && v(1,i)+step>0)
           l(i) = 1;
        end
    end 
    
    if(~isempty(find(l,1)))
        
        % If there is two variables that contains 0, we pick one randomly.
        s = find(l,1); 
    
        % Here, we need to detect with monomials containes the variables. 
        ind = Expo(:,s)>0;
        tmp = [tmp tmp(:,1)]; 
        tmp(ind,:) = 0;  
    end
    
    i1 = find(Coeff>=0); i2 = find(Coeff<0);
    t1 = max(tmp(i1,:),[],2); t2 = min(tmp(i2,:),[],2);

    K = t1'*Coeff(i1)+t2'*Coeff(i2);

     if(strcmp(type,'min'))
         K = -K;
     end

     % K = sqrt(K)*norm(step*sqrt(dim/4),2); 
     
end
            