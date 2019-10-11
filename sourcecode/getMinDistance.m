
function minDist=getMinDistance(X,Center,dim_epi)
       
       minDist=dim_epi;
       if isempty(Center)
           minDist=dim_epi;
       else
          L=length(Center(:,1));
          for i=1:L
              c=dim_epi - length(intersect(Center(i,:),X));
%               c = dim_epi;
%               for k = 1:dim_epi
%                   for s = 1:dim_epi
%                       if X(k) == Center(i,s)
%                           c = c + 1;
%                       end
%                   end
%               end
              if c<minDist
                  minDist=c;
              end              
          end
       end     
 end