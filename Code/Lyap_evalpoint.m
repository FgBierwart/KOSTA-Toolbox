
% This function computes the lyapunov function on the grid defined by c
% (see main.m for the construction of c). 

function [V,dot_V] = Lyap_evalpoint(basis,Vec,indx,c,s,approx)

n = s.n; dim = size(c,1); F = s.field; 
fct = ones(n,length(c));

if(strcmp(basis,'monomials'))
    combs = s.Expo;
    for i = 1:dim
        fct = fct.*(c(i,:).^combs(:,i));  
    end
else
    if(approx.basis~=0)
        temp = cell(1,dim);
        h = 1e-10; 
        for k = 1:dim
            tmp = c'; 
            tmp(:,k) = tmp(:,k)+h;
            temp{k} = tmp;
        end
    end
    M = cell(dim,1);
    for i = 1:s.n
        if(approx.basis==0)
            temp = s.grad_f{i}(c'); 
        end
        fct(i,:) = s.f{i}(c'); 
        for k = 1:dim
            if(approx.basis==0)
                M{k} = [M{k};temp(:,k)'];
            else
                tmp = (s.f{i}(temp{k})'-fct(i,:))./h;
                M{k} = [M{k};tmp];
            end
        end
    end
end

eval_F = F(c');
sol = (fct')*Vec(:,indx); coeff = ones(1,dim); 
V = sum((abs(sol).^2)*coeff',2);

dot_V = zeros(length(c),1);

for k = 1:dim
    if(strcmp(basis,'monomials'))
        Mat = zeros(n,length(c));
        temp = combs;
        temp(:,k) = temp(:,k)-1;
        
        h_n = find(temp(:,k)<0); 
        h_z = find(temp(:,k)>=0 & sum(temp,2)==0); 
        h_p = find(temp(:,k)>=0 & sum(temp,2)~=0);  
        
        t = temp(h_p,:);
        [~,b] = ismember(t,combs,'rows');
        
        Mat(h_p,:) = (t(:,k)+1).*fct(b,:); 
        Mat(h_z,:) = ones(length(h_z),length(c));
        Mat(h_n,:) = zeros(length(h_n),length(c));
    else
        Mat = M{k}; 
    end
    
    T = Mat'*Vec(:,indx);
    dot_V = dot_V + ...
             sum(real(T).*real(sol) + imag(T).*imag(sol),2).*eval_F(:,k);
end

end
