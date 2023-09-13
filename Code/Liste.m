
% Given a grid 'liste' (as final_grid in AGM.m for instance), this function 
% uptade the grid by dividing each cells in (div)^n parts. 'liste' contains 
% the final result. 

function [liste] = Liste(liste,div)
    
    dim = size(liste(:,1:(end-1)),2);

    CC = cell(1,dim);
    [CC{:}] = ndgrid(0:(div-1));
    combs = [];
    for i = 1:length(CC)
        combs = [combs CC{i}(:)];
    end

    idx = size(unique(liste(:,1:(end-1)),'rows'),1);
    sol = repmat(combs,idx,1);

    delta_x = liste(:,end);

    Temp = repelem(liste(:,1:(end-1)),div^dim,1) + ...
        sol.*(repelem(delta_x,div^dim,1)/div);

    clear liste
    liste = zeros(size(Temp,1),size(Temp,2)+1);
    liste(:,1:(end-1)) = Temp;
    liste(:,end) = repelem(delta_x,div^dim,1)/div;

end