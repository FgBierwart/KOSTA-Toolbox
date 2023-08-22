
% -------------------------------------------------------------------------
% this function provides a first reference level for the bisection method
% in the extraction of the region of attraction (see Lyap_certificate.m)
% -------------------------------------------------------------------------

function gam2 = ref_level_set(V,nablaV,f,err,Vmin,Vmax,p,dim,approx)

x = mpvar('x',1,dim);
Z2 = monomials(x,0:p);
options.solver = 'mosek';  
tol = 0.5; 
r = sqrt(2); 
d = sort(linspace(Vmin,Vmax,10),'descend');
k = 1; 

combs = Exp_mon(dim,dim-1); [a,~] = find(combs>1); combs(a,:)=[];
combs = [combs;ones(1,dim)]; combs(combs==0)=-1; combs = combs'; 

while(true)
    
    gam2 = d(k);
    gam1 = gam2;
    
    if(approx.flag~=0)

        flag = zeros(1,2^dim); 
    
        for i = 1:(2^dim)
            prog = sosprogram(x);
            for j = 1:(dim+3)
                [prog,s{j}] = sossosvar(prog,Z2);
            end
            l = 1e-08*(x*x'); 
            q = (nablaV'*f)+(nablaV.*combs(:,i))'*err+l;
            %prog = sosineq(prog,-q-s{end}*(r-sum(x.^2))+s{1}*(V-gam2)-s{2}*(V-gam1)-[s{dim+1:end-1}]*(nablaV.*combs(:,i)));
            prog = sosineq(prog,-q-s{end}*(r-sum(x.^2))+s{1}*(V-gam2)-s{2}*(V-gam1));
            [~,info] = sossolve(prog,options);
            flag(i) = info.feasratio;
        end

    else
        prog = sosprogram(x);
        for j = 1:(dim+1)
            [prog,s{j}] = sossosvar(prog,Z2);
        end
        l = 1e-06*(x*x'); 
        q = (nablaV'*f+l);
        prog = sosineq(prog,-q-s{end}*(r-sum(x.^2))+s{1}*(V-gam2)-s{2}*(V-gam1));
        [~,info] = sossolve(prog,options);
        flag = info.feasratio;
    end

    if(flag>=tol)
        break; 
    end

    k = k+1;
    disp(gam2);

end

end