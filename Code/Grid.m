
% -------------------------------------------------------------------------
% This function generates a matrix of data points in 'dim' dimension where 
% each row is associated to a coordinate. 
% -------------------------------------------------------------------------

function c = Grid(a,delta,n,dim)

    xr = cell(dim,1); 
    
    for i = 1:dim
        xr{dim-i+1} = linspace(a(i),a(i)+delta,n); 
    end
            
    dim = length(xr);
    x = cell(1,dim);
    [x{:}] = ndgrid(xr{:});
    
    c = [];
    
    for i = dim:-1:1
        X = x{i};
        c = [c;X(:)'];
    end

end
