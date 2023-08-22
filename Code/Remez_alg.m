function [alpha,h] = Remez_alg(f,n,z,a)

    syms x1 x2;
    fct = matlabFunction(f,'Vars',{[x1;x2]});    
    Z = matlabFunction(z,'Vars',{[x1;x2]});
    
    % Construction of the domain (r sample points) 

    % 1. Random nodes
    
    % r = n+1; 
    % c = (a-(-a))*rand(2,r)-a; 

    % 2. Chebyshev nodes
    
    A = [-a,-a]; B = [a, a]; 
    c = computeChebyshevNodes2D(A,B,50)'; 
    r = size(c,2);
    
    % Contruction of Gram matrix 
    
    A = zeros(r,n); 

    for i = 1:size(A,2)
        tmp = matlabFunction(z(i),'Vars',{x1,x2}); 
        A(:,i) = tmp(c(1,:),c(2,:));
    end
    
    F = fct(c);
    c1 = [F -F]; 
    tmp = ones(1,r); 
    C1 = [A' -A';tmp tmp];
    k = 0; 
    h = 100; 
    err = matlabFunction(sym(0),'Vars',{x1,x2});  

    while(true)

        c_tmp = (a-(-a))*rand(2,10000)-a; 
    
        if(k>=1)
    
            err = matlabFunction(abs(f-value(alpha)'*z));
            tmp = sort(err(c_tmp(1,:),c_tmp(2,:)),'descend');
            [~,id] = intersect(err(c_tmp(1,:),c_tmp(2,:)),tmp(1:1));
            vec = c_tmp(:,id); 

            d = zeros(size(C1,1)-1,2*size(vec,2)); 
            d(:,1:2:end) = Z(vec); 
            d(:,2:2:end) = -Z(vec);
            C1 = [C1(1:end-1,:) d;C1(end,:) ones(1,size(d,2))];

            F = zeros(1,2*size(vec,2)); 
            F(:,1:2:end) = fct(vec); 
            F(:,2:2:end) = -fct(vec);
            c1 = [c1 F]; 
            
        end
        
        if(abs(max(err(c_tmp(1,:),c_tmp(2,:)))-value(h))<1e-04)
            break;
        end
                            
        alpha = sdpvar(n,1);
        sdpvar h; 
        con = [C1'*[alpha;h]-c1'>=0]; 
        diagnostics = solvesdp(con,h,sdpsettings('verbose',1)); 
        k = k+1; 
        
    end

end

function xi = computeChebyshevNodes2D(a, b, N)
    t = cos(pi * ((0:N-1)' + 0.5) / N); % Chebyshev nodes in one dimension
    
    xi1 = (b(1) - a(1)) / 2 * t + (b(1) + a(1)) / 2; % Chebyshev nodes in first dimension
    xi2 = (b(2) - a(2)) / 2 * t + (b(2) + a(2)) / 2; % Chebyshev nodes in second dimension
    
    [xi1, xi2] = meshgrid(xi1, xi2); % Create grid of 2D Chebyshev nodes
    xi = [xi1(:), xi2(:)]; % Combine the nodes into a single matrix
end