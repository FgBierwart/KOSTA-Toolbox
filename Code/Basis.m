function [s,approx] = Basis(dim,choice,approx)
    
    s = struct;    
    x = sym('x',[1 dim]); 

    if(strcmp(choice,'monomials'))

        % Monomials 
        % ---------
        
        order = 5; 
        combs = Exp_mon(dim,order);
        combs(1,:) = [];
        s.n = size(combs,1);
        s.Expo = combs;
        approx.basis = 0; 

    else

        % Gaussian functions
        % ------------------
        
        n = 9; a = 1; gamma = 0.1; 
        s.f = cell(1,n);  
        ci = Grid(-a*ones(dim,1),2*a,nthroot(n,dim),dim); 
        s.grad_f = s.f; 

        for k = 1:n
            temp = exp(-(gamma^2)*sum((x-ci(:,k)').^2,2));
            s.fsym(k) = temp; 
            s.grad_f{k} = matlabFunction(transpose(gradient(temp)),'Vars',{x});
            s.f{k} = @(x) exp(-(gamma^2)*sum((x-ci(:,k)').^2,2)); 
        end

        s.n = n;  
        approx.basis = 0;
        
    % Complete for other basis functions  
    % ... 

    end

end