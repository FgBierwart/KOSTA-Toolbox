
% =========================================================================
% INPUT: 
% -----
% - basis : string caracter that contains the basis.
% 
% - s : structure that contains all information about the basis of
%       interest. E.g : 
%
%       * For monomials: s contains the order of monomials (Expo)
%       * For gaussian: s contains each basis functions encoded in s.f. 
%
% - Vec : Matrix of approximate eigenvectors of the Koopman operator
% 
% - indx : column of Vec that indicates the eigenvectors with closest 
%          eigenvalues of the one of jacobian matrix.
%
% - F_exp : Cell encoding the polynomial vector field or its approximation.
%           Each cell is a matrix associated to a componenent of the
%           vector. Each rows contains monomials and associated
%           coeffificient. E.g, F = [x(2) -2*x(1)+(1/3)*x(1)^3-x(2)];
% 
%           F_exp{1} = 0 1 1
%           F_exp{2} = |0 1 -1 |
%                      |1 0 -2 |
%                      |3 0 1/3|
%
%
% - approx : structure that contains: approx.flag which is 0
%            for polynomial field or 1 for nonpolynomial vector fields, 
%            approx.err encodes the polynomial error (for non-polynomial 
%            vector field) in a same way as F_exp for each componenent of 
%            the vector field (E.g with taylor: err = c(i)||x||^d.) For 
%            polynomial vector fields, error is 0 for each component. It
%            also contains approx.basis which is 1 if the gradient of the
%            basis functions are approximated and 0 if there are given by
%            the user.    
%             
% 
% - p : degree of multipliers in SOS. 
%
% OUTPUT: 
% ------
% 
% This function return the basin of attraction approximation delimited by
% the set {gam_1 < V < gam_2} and computed using SOS validation method.
% Here V is polynomial variable representing the Lyapunov function. 
%
% =========================================================================

function [gam1,gam2,Vmin,V] = Lyap_certificate(basis,s,Vec,indx,F_exp,approx,p)

% =========================================================================
% Construction of V
% ========================================================================= 

dim = length(indx); 
choice = 'minimax'; a = 1; 

% Computation of the polynomial lyapunov function (exact with monomials or
% approximation using a minmax approximation for gaussian functions). 

if(~strcmp(basis,'monomials')) % minimax approx of the lyapunov function 
    
    order = 12;    
    Coeff = [real(Vec(:,indx)) imag(Vec(:,indx))]; 

    if(strcmp(choice,'minimax'))

        % First possibility, we construct the Lyapunov approximation based
        % on the different evaluation functions 

        if(approx.basis~=0)

            Vg = struct;
            Vg.basis = s.f;
            Vg.coeff = Coeff;
            flag = 1; 

        else

        % Or we compute the symbolic expression and we obtain the remez
        % version of the symbolic expression directly. (default) 
        
        g = s.fsym; 
        Vg = sum(sum(Coeff.*transpose(g),1).^2);
        flag = 0; 

        end

    elseif(strcmp(choice,'Taylor'))
         g = s.fsym; 
         Vg.expr = sum(sum(Coeff.*transpose(g),1).^2); % symbolic
         flag = 0; 
    end
        
    [Coeff,Expo,~] = Pol_approx(Vg,a,choice,order,dim,flag);
    
    % tmp = find(abs(Coeff)<1e-03); 
    % Coeff(tmp) = []; Expo(tmp,:) = [];  

else
    [Coeff,Expo] = Lyapunov_pol(s.Expo,Vec,indx);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Here we are computing the highest level set of V that is contained in the
% [-1,1]^2. This will be the starting point of the bisection ensuring that
% no levelset is picked outside the domain of interest. 

% 2D 
N = 1000; t = Hypercube_data(dim,N,-0.98*ones(dim,1),2);

% 3D
% N = 625; t = Hypercube_data(dim,N,-ones(dim,1),2);

P = eval_pol(Coeff,Expo,t);
Vmax = min(P); 
disp(Vmax);
% -------------------------------------------------------------------------

% =========================================================================
% 2. Application of SOS
% =========================================================================

% constructiox = mpvar('x',1,dim); n of the lyapunov candidate
x = mpvar('x',1,dim); 
V = (Coeff(1))*prod(x.^Expo(1,:),2);

for i = 2:length(Coeff)
    V = V + Coeff(i)*prod(x.^Expo(i,:),2);
end

nablaV = []; 
for i = 1:dim
    nablaV = [nablaV;diff(V,x(i))];
end

f = mpvar('p',[dim,1]); 

% vector field 
for i = 1:dim
    f(i) = (F_exp{i}(1,end))*prod(x.^F_exp{i}(1,1:end-1),2);
    for j = 2:size(F_exp{i},1)
        f(i) = f(i) + (F_exp{i}(j,end))*prod(x.^F_exp{i}(j,1:end-1),2); 
    end
end

% error made on nonpolynomial vector field 
if approx.flag~=0
    err = mpvar('e',[dim,1]); 
    for i = 1:dim
        err(i) = (approx.err.Coeff{i}(1))*prod(x.^approx.err.Expo{i}(1,:),2);
        for j = 2:length(approx.err.Coeff{i})
            err(i) = err(i) + approx.err.Coeff{i}(j)*prod(x.^approx.err.Expo{i}(j,:),2);
        end
    end
else
    err = mpvar('e',[dim,1]);
    for i = 1:dim
        err(i) = 0;
    end
end

% 4. If we approximate V, we are going to compute the minimum value so that
% during the bisection, we can not go below this values since it does not
% make sens. 
% -------------------------------------------------------------------------

if(strcmp(basis,'gaussian'))
    options.solver = 'mosek'; 
    c = mpvar('c',1,2*dim); k = 1;  
    for i = 1:dim
        for j = 0:1
            c(k) = 1+(-1)^j*x(i); 
            k = k+1; 
        end
    end
    [Vmin,~,~] = findbound(V,c,[],order,options);
    V0 = eval_pol(Coeff,Expo,zeros(1,dim));
    while(abs(Vmin-V0)<1e-06)
        Vmin=1.1*Vmin;
    end
else
    Vmin = 0;
end

% bisection to extract the BOA 
Vref = ref_level_set(V,nablaV,f,err,Vmin,Vmax,p,dim,approx);
gam2 = bisection_top(V,nablaV,f,err,Vref,Vmax,p,dim,approx);
gam1 = bisection_bot(V,nablaV,f,err,Vref,Vmin,p,dim,approx);

end