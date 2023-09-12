
% =========================================================================
% INPUT: 
% -----
%
% - dim : number of variables in the function we want to approximate
% 
% - f : strutcure that encode the function we want to approximate.
% 
%        * if flag = 0, f is a symbolic function. 
%        * if flag = 1, f contains f.basis (cell of handle functions) and 
%          f.coeff (Coefficient). For instance, 
%
%          g(x) = exp(x(1))-2*cos(x(2)) is encoded by 
%
%               o f.basis = {@(x) exp(x(1)),@(x) cos{x(2)}},
%               o f.coeff = [1,-2].
% 
% - n : number of basis function 
%
% - z : vector of symbolic basis function used to approximate the function
%       encoded by f. For instance, monomials up to order 5. 
% 
% - a : encdes the domain over which the function is approximated. The
%       function is approximated over X = [-a,a]^dim.
%
%
% OUTPUT: 
% ------
% 
% - alpha : Coefficient associted to variable in z resulting for the 
%           minimax approximation of f over X using the basis z. 
% 
% - h : error of approximation. If f_app is the approximation of f, h is
%       such that |f-f_app| < h. (Solution of relaxed LP program. See the 
%       paper for more detailed)
%
% =========================================================================

function [alpha,h] = Remez_alg(f,n,z,a,flag,dim)
    
    x = sym('x',[1 dim]);
    Z = matlabFunction(z,'Vars',{x});
    
    if(flag==0)
        f = matlabFunction(f,'Vars',{x}); 
    end
    
    % Construction of the domain (Chebyshev nodes)
    
    r = 50; 
    A = -a*ones(1,dim); B = -A; 
    c = computeChebyshevNodesND(A,B,r,dim); 
    r = size(c,2);
    
    % Contruction of Gram matrix 
    
    A = zeros(r,n); 

    for i = 1:size(A,2)
        tmp = matlabFunction(z(i),'Vars',{x}); 
        A(:,i) = tmp(c');
    end
    
    if(flag~=0)
        for i = 1:size(f.coeff,1)
            temp(i,:) = f.basis{i}(c'); 
        end
        F = sum((f.coeff'*temp).^2,1); 
    else
        F = f(c')';
    end
    
    c1 = [F -F]; 
    tmp = ones(1,r); 
    C1 = [A' -A';tmp tmp];
    k = 0; 
    h = 100; 
    err = matlabFunction(sym(0),'Vars',{x});
    err = err(c'); 
    err2 = 10; 

    % This loop serves to update the domain of chebysehv nodes by adding
    % points over a random domain such that error of approximation is the 
    % largest one. The algorithm stop when the difference two succesive 
    % error is smaller than 1e-04.    

    while(true)

        c_tmp = (a-(-a))*rand(dim,10000)-a;
    
        if(k>=1)

            if(flag~=0)
                clear temp; 
                for i = 1:size(f.coeff,1)
                    temp(i,:) = f.basis{i}(c_tmp'); 
                end
                R = sum((f.coeff'*temp).^2,1)';
            else
                R = f(c_tmp');
            end
            
            app = matlabFunction(value(alpha)'*z,'Vars',{x});
            err2 = abs(R-app(c_tmp'));
            tmp = sort(err2,'descend');
            [~,id] = intersect(err2,tmp(1:1));
            vec = c_tmp(:,id); 

            d = zeros(size(C1,1)-1,2*size(vec,2)); 
            d(:,1:2:end) = Z(vec'); 
            d(:,2:2:end) = -Z(vec');
            C1 = [C1(1:end-1,:) d;C1(end,:) ones(1,size(d,2))];

            F = zeros(1,2*size(vec,2));

            if(flag~=0)
                clear temp; 
                for i = 1:size(f.coeff,1)
                    temp(i,:) = f.basis{i}(vec'); 
                end
                F(:,1:2:end) = sum((f.coeff'*temp).^2,1);
                F(:,2:2:end) = -F(:,1:2:end); 
            else
                F(:,1:2:end) = f(vec'); 
                F(:,2:2:end) = -f(vec');
            end
            
            c1 = [c1 F]; 
            
        end
        
        if(abs(max(err)-max(err2))<1e-04)
            break;
        else
            err = err2; 
        end
                            
        alpha = sdpvar(n,1);
        sdpvar h; 
        con = [C1'*[alpha;h]-c1'>=0]; 
        diagnostics = solvesdp(con,h,sdpsettings('verbose',0)); 
        k = k+1;
        
    end

end

% This function computes N^dim Cehbyshev nodes in dim dimension over 
% X = [a(1),b(1)] x ... x [a(n),b(n)]. 

function xi = computeChebyshevNodesND(a,b,N,dim)

    t = cos(pi * ((0:N-1)' + 0.5) / N); % Chebyshev nodes in one dimension

    for i = 1:dim
        xr{i} = (b(i) - a(i)) / 2 * t + (b(i) + a(i)) / 2; % Chebyshev nodes in ith dimension
    end

    x = cell(1,dim);
    [x{:}] = ndgrid(xr{:});
    
    xi = [];
    
    for i = dim:-1:1
        X = x{i};
        xi = [xi;X(:)'];
    end
    
end