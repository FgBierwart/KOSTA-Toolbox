
% According to the polynomial describe by Expo and Coeff, we are computing 
% the gradient of the polynomial. Expo_g is a cell where each cell is
% associated to one componenent of the gradient. 

function [Expo_g,Coeff_g] = Gradient_pol_mon(Expo,Coeff)

dim = size(Expo,2);

    for i = 1:dim
       t = Expo;
       Expo_g{i} = Expo; 
       Coeff_g{i} = Coeff;
       Expo_g{i}(:,i) = Expo_g{i}(:,i) - 1; 
       temp = find(Expo_g{i}(:,i)==-1);
       t(temp,:) = [];
       Expo_g{i}(temp,:) = []; 
       Coeff_g{i}(temp,:) = [];
       Coeff_g{i} = Coeff_g{i}.*t(:,i);

       H = unique(Expo_g{i},'rows');
       D = zeros(size(H,1),size(Coeff_g{i},2));
    
        for k = 1:size(H,1)
            id = ismember(Expo_g{i},H(k,:),'rows');
            D(k,:) = sum(Coeff_g{i}(id,:),1);   
        end

        Expo_g{i} = H;
        Coeff_g{i} = D;

    end    

end