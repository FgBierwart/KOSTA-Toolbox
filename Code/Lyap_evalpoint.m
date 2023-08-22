
% This function computes the lyapunov function on the grid defined by c
% (see main.m for the construction of c). This fnction is implemented now
% for gaussian and monomials. 

function [V,dot_V] = Lyap_evalpoint(basis,Vec,indx,c,s)

n = size(Vec,1); dim = size(c,1); F = s.c; 
fct = ones(n,length(c));

if(strcmp(basis,'monomials'))
    combs = Exp_mon(dim,s.a); combs(1,:)=[]; 
    for i = 1:dim
        fct = fct.*(c(i,:).^combs(:,i));  
    end
elseif(strcmp(basis,'gaussian'))
    gamma = s.a; 
    don = s.d'; 
    for i = 1:n
         t = exp((-gamma^2)*sum((c-don(i,:)').^2,1));
         fct(i,:) = t;
    end
end

eval_F = F(c');
sol = (fct')*Vec(:,indx); coeff = ones(1,dim); 
V = sum((abs(sol).^2)*coeff',2);

Mat = zeros(n,length(c));
dot_V = zeros(length(c),1);

for k = 1:dim
    if(strcmp(basis,'monomials'))
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
    elseif(strcmp(basis,'gaussian'))
        Mat = (-2*gamma^2).*fct.*(c(k,:)-don(:,k));
    end
    T = Mat'*Vec(:,indx);
    dot_V = dot_V + ...
             sum(real(T).*real(sol) + imag(T).*imag(sol),2).*eval_F(:,k);
end

end
