
% This function return 0 if the level curve does not cross the cells and 1 
% otherwise. 

function flag = CrossingCurve(G,q,V,V_star,C,type,K,P) 

flag = zeros(length(q),1); 

if(isempty(q))

    flag = 0; 

else
    
    if(sum(G(q,end))==length(q))
        
        count_tc = [];
        
        for i = 1:length(q)
            
            flag(i) = CrossingCurve_pct(G,V,q(i),V_star,C,type,K(q(i)),P); 
            
            % The criterion can not ensure we are not crossing the cell
            if(flag(i)~=0)
                % so we compute the neighbor of the cell involve
                tmp = Neighbour(G,q(i));
                if(strcmp(type,'top'))
                    w = min(V(tmp,:),[],2)>V_star;
                elseif(strcmp(type,'bot'))
                    w = max(V(tmp,:),[],2)<V_star;
                end
                tmp = tmp(w,:);
                count_tc = unique([count_tc;tmp]); 
            end           

        end     

        [~,b] = intersect(count_tc,q);
        count_tc(b) = []; 

        if(isempty(find(flag==1,1)))

            clear flag; 
            flag = 0; 

        else 

            if(sum(G(count_tc,end))==length(count_tc))

                flag = zeros(length(count_tc),1);

                for i = 1:length(count_tc)
                    flag(i) = CrossingCurve_pct(G,V,count_tc(i),V_star,C,type,K(count_tc(i)),P);
                end

                if(isempty(find(flag==1,1)))
                    clear flag; 
                    flag = 0;
                else
                    clear flag; 
                    flag = 1; 
                end
                
            else
                flag=1;
           end

        end

    else
        flag = 1;
    end

end

end

