
function flag = CrossingCurve_pct(G,V_grid,indx,V_star,C,type,K,P)

Expo = C{1}; Coeff = C{2}; 
% E_nrm_grad = C{3}; C_nrm_grad = C{4};
v = C{5}; 

N = 700; t = Hypercube_data(size(Expo,2),N,G(indx,:)',G(indx,end-1));
K = sqrt(K)/(2*N); 

P = P(t); 

flag = 1; 

if(strcmp(type,'top'))
    Kmin = Bound_Grad_Pol(Expo,-Coeff,v{indx},G(indx,end-1),'min'); 
    vec = min(P)-K; 
    if(min(V_grid(indx,:))>V_star && (vec>V_star || Kmin>V_star))
        flag = 0;
    end
elseif(strcmp(type,'bot'))
    vec = max(P)+K; 
    Kmax = Bound_Grad_Pol(Expo,Coeff,v{indx},G(indx,end-1),'max'); 
    if(max(V_grid(indx,:))<V_star && (vec<V_star || Kmax<V_star))
        flag = 0;  
    end
end
