function [P_value, G] =Chi_square2(snp_com,state)
%compute the chi-score
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end


[xrow,xcol] = size(snp_com);

[Data,idx,cid]=unique(snp_com,'rows');
[lrow,~]=size(Data);



%% 

F=zeros(3,lrow+1);  %% ��һ�д�� case ���� �ڶ��д�� control ���� �����д�� ����
E=zeros(2,lrow);  %% ������

for i=1:xrow   %% ͳ��ÿ����������ϳ��ֵĴ���
   
       if state(i)==0 %% case        
            F(1,cid(i))=F(1,cid(i))+1;        
       else 
            F(2,cid(i))=F(2,cid(i))+1;  
       end  
   
end

F(3,1:lrow)=sum(F(1:2,1:lrow),1);
F(:,lrow+1)=sum(F(:,1:lrow),2);


G=0;
Degree=(2-1)*(lrow-1);%(2-1)*(lrow-1)

for i=1:2
    for j=1:lrow
        O=F(i,j);
        E=( F(i,lrow+1) * F(3,j) )/xrow ;
       
     
          if E>0
              G=G+(abs(O-E))^2/E;
           end
     
      
    end    
end
%Degree
G=2*G;

P_value=1-chi2cdf(G,Degree);


