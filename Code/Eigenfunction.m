
% =========================================================================
% INPUT: 
% -----
% - s : structure that contains all information about the basis of
%       interest. E.g : 
%
%       * For monomials : s contains the order of monomials
%       * For gaussian, s contains the value of gamma, the domain of the 
%         centers and the centers themselves.
%
% - basis : string caracter that contains the basis.
%
% - eig_ex : exact eigenvalues of the vector fields chosen to compute the
%            lyapunov function (either polynomial or its polynomial 
%            approximation. See Readme (flowchart) for more details.
%
% - truncation : indicates the type of projection considered. truncation =
%                1 for a truncation or 0 for orthogonal projection. 
%
%
% OUTPUT: 
% ------
% - Vec : Matrix of approximate eigenvectors of the Koopman operator
% 
% - indx : column of Vec that indicates the eigenvectors with closest 
%          eigenvalues of the one of jacobian matrix.
% 
% - Val_p : Diagonal matrix containing approximate eiganvalues.
% 
% - L : approximation of the infinitesimal generator.
% =========================================================================

function [Vec,indx,Val_p,L] = Eigenfunction(s,basis,eig_ex,truncation) 

% Computation of infinitesimal generator matrix approximation : 

if(truncation==1)
    
    % =====================================================================
    % Truncation (monomials)
    % ===================================================================== 
    
    dim = length(s.b); 
    order = s.a;
    F = s.b; 
    
    combs = Exp_mon(dim,order);
    combs(1,:) = [];
    n = size(combs,1);
    
    L = zeros(n,n);
    
    for i = 1:n
    
        Expo = cell(dim,1);
        t = combs(i,:);
        
        for k = 1:dim
            tmp = zeros(1,dim);
            tmp(k) = 1;
            Expo{k} = [F{k}(:,1:dim)+(t-tmp) (t(k)).*F{k}(:,end)];
        end
        
        A = cat(1,Expo{:});
        
        H = unique(A(:,1:end-1),'rows');
        H(:,dim+1) = zeros(size(H,1),1);
        
        for k = 1:size(H,1)
            id = ismember(A(:,1:dim),H(k,1:dim),'rows');
            H(k,end) = sum(A(id,end));   
        end
        
        [~,b] = ismember(H(:,1:end-1),combs,'rows');
        H(b==0,:) = [];
        b(b==0)=[];
        L(b,i) = H(:,end);
    
    end

    L = L'; 

else

    % =====================================================================
    % L2 Projection (Gaussian or monomials) 
    % =====================================================================
    
    if(strcmp(basis,'gaussian'))

        gamma = s.a; n = s.b; F = s.c; don = s.d'; dim = size(don,2);
        a = 0.1; m = 64; xx = Grid(-a*ones(dim,1),2*a,nthroot(m,dim),dim); 
        eval_F = F(xx');    
        Psi_y = zeros(n,length(xx));
        Psi_x = zeros(n,length(xx));
        
        for i = 1:n
            t = exp((-gamma^2)*sum((xx-don(i,:)').^2,1));
            Psi_x(i,:) = t;
            for k = 1:dim
                Psi_y(i,:) = Psi_y(i,:) + Psi_x(i,:).*eval_F(:,k)'.*(-2*gamma^2).*(xx(k,:)-don(i,k)); 
            end    
        end

    elseif(strcmp(basis,'monomials'))

        F = s.c; order = s.a; dim = length(s.b); 
        combs = Exp_mon(dim,order); combs(1,:) = []; n = size(combs,1);
        a = 0.2; m = 50; xx = Grid(-a*ones(1,dim),2*a,m,dim); 
        eval_F = F(xx');

        Psi_y = zeros(n,length(xx));
        Psi_x = zeros(n,length(xx));

        for k = 1:dim
            Expo_g{k} = combs; 
            Coeff_g{k} = Expo_g{k}(:,k);
            Expo_g{k}(:,k) = Expo_g{k}(:,k)-1; 
            tmp = Expo_g{k}(:,k)==-1; 
            Coeff_g{k}(tmp) = 0; 
            Expo_g{k}(tmp,k) = 0; 
        end
        
        for i = 1:n
            t = prod(xx.^(combs(i,:)'));
            Psi_x(i,:) = t;
            for k = 1:dim
                Psi_y(i,:) = Psi_y(i,:) + Coeff_g{k}(i)*prod(xx.^(Expo_g{k}(i,:)')).*eval_F(:,k)'; 
            end    
        end
    end 
    
    L = (Psi_y*Psi_x')*(pinv(Psi_x*Psi_x'));
    %L = Psi_y*pinv(Psi_x); 

end

% =========================================================================
% Computation of eigenfunctions
% =========================================================================

[Vec,Val_p] = eig(L');
indx = zeros(dim,1); 

% here, we are extrating the eigenvalues that are closest to the jacobian. 

for i = 1:length(eig_ex)

    D = diag(Val_p);
    temp = abs(eig_ex(i)-D);
    d = sort(temp);   
    lowest = find(temp==d(1));
    D(lowest) = Inf;
    indx(i) = lowest;

end

% here, we remove almost zeros coefficients for stability numerical purpose
% (small terms in the Lyapunov candidate and since it a candidate we do not
% lose generality).  

% if strcmp(basis,'monomials') 
%     S = abs(Vec(:,indx));
%     for i = 1:size(S,2)
%         tmp = find(S(:,i)<1e-04); 
%         Vec(tmp,indx(i)) = 0; 
%     end
% end

end
