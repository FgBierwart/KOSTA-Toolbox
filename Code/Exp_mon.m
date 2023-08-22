
% -------------------------------------------------------------------------
% This function generates exponants of monomial up to order 'order'. For
% instance, if order = 2: 
% 
% combs : 0 0
%         0 1
%         1 1
%         2 0
%         0 2
% -------------------------------------------------------------------------

function combs = Exp_mon(dim,order)

    Exp = 0:order;
    C = cell(1,dim);
    
    for i = 1:dim
        C{i} = Exp;
    end
    
    combs = cell(1, dim);
    [combs{:,:}] = ndgrid(C{1:dim});
    combs = reshape(cat(dim+1,combs{:}), [], dim);
    temp = find(sum(combs,2) > order);
    
    
    for i=1:length(temp)
        combs(temp(i),:) = [];
        temp((i+1):length(temp),:) = temp((i+1):length(temp),:)-1;
    end
    
    combs = combs(:,end:-1:1);

end

