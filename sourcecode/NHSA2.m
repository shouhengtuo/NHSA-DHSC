function [Candidate,canSize,NC,totaltime,flag] = NHSA2(data,dim_epi,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX)
%input--------------------------------------------------------------------
% data-----------------input dataset
% epi_dim--------------the epistasis dimension
% HMS--------------the size of harmony memory(HM)



%%-------------------------------------------------------------------------
% initial arguments
HMCR=0.96;
PAR=0.35;
n=size(data,2);
State=data(:,n);


Candidate1=ones(CandidateSize,dim_epi+3); %% 存放候选解
canSize1=0;
Candidate2=ones(CandidateSize,dim_epi+3); %% 存放候选解
canSize2=0;
Candidate3=ones(CandidateSize,dim_epi+3); %% 存放候选解
canSize3=0;

%% ---------------------------------------------------------------

SNPs=n-1;  %% 总SNP个数



%% 初始搜索集合
for i = 1: SNPs 
    S_pvalue(i) = Gtest_score(data(:,i), State);
end

[S_P, SINDEX] = sort(S_pvalue);

% fprintf('%d ',SINDEX(1:HMS))
% fprintf('\n');




%%

EliteSize=fix(HMS/5);  %% HMS的偶数倍
Elite1=[]; %% 存放精英解集。
Efit1=[];

Elite2=[]; %% 存放精英解集。
Efit2=[];

Elite3 = []; %% 存放精英解集。
Efit3 = [];


Center=zeros(1,dim_epi+1); %% 存放精英中心集
% maxIterForLocalSearch=min(fix(max_iter/5),3000);  %% 每次循环最大值，超过时，记录中心和精英，清空解集


%% 初始化
flag = -1;
X=zeros(HMS,dim_epi);
snp=[];
for i=1:HMS
%     RR = randperm(HMS);
%      snp(1:dim_epi) = SINDEX(RR(1:dim_epi));
%%    
%     snp(1)=SINDEX(i);
%     for j=2:dim_epi
%       snp(j)=ceil(rand*SNPs);  
%     
%       while ismember(snp(j),snp(1:j-1)) 
%          snp(j)=ceil(rand*SNPs);        
%       end
%     end
%%
     snp(1)=ceil(rand*SNPs);
    for j=2:dim_epi
      snp(j)=ceil(rand*SNPs);  
    
      while ismember(snp(j),snp(1:j-1)) 
         snp(j)=ceil(rand*SNPs);        
      end
    end
%%

   temp=snp;
    snp=sort(snp);
    while ismember(snp,X,'rows')
        j=ceil(rand*dim_epi);
        snp(j)=ceil(rand*SNPs); 
        temp=snp;
        snp=sort(snp);
    end
      
    X(i,:)=snp;   %% X中存放有序的解
    if snp == CX
        flag = 111;
    end
    HM(i,:)=temp;  %% HM中相应存放无序解
    [Fit(i,1),Fit(i,2),Fit(i,3)] = ScoreFunctions(data(:,X(i,:)),State);    
    snp=[];
    
%     fprintf('[%d %d]  ', X(i,1),X(i,2));
end
fprintf('\n');
 



X2=X;
HM2=HM;
Fit2=Fit;

X3 = X;
HM3 = HM;
Fit3 = Fit;


NC=0;
 LT=0;
%%-------------------------------------------------------------------------
tic;
while NC <= max_iter      
    flag = 0;
    [best1,idbest1] = min(Fit(:,1));
    [best2,idbest2] = min(Fit2(:,2));
    
    [best3,idbest3] = min(Fit3(:,3));
   
    Xbest1 = HM(idbest1,:);
    Xbest2 = HM2(idbest2,:);
    Xbest3 = HM3(idbest3,:);
    IDBEST = [idbest1, idbest2, idbest3];
     i=1;
     while i<=dim_epi 
         if rand<HMCR
                a=ceil(rand*HMS*3);
%                 r = ceil(rand*3);
%                 a = IDBEST(r);
                if a>2*HMS
                    Xnew(i)=HM(a - 2*HMS,i); 
                    if rand<PAR                    
                       Xnew(i)=Xnew(i)+((-1)^Xnew(i))*abs(Xbest1(i)-HM(ceil(rand*HMS),i));
                        Xnew(i)=max(min(Xnew(i),SNPs),1);
                    end   
                elseif a>HMS
                     Xnew(i)=HM2(a-HMS,i);
                     if rand<PAR                     
                        Xnew(i)=Xnew(i)+((-1)^Xnew(i))*abs(Xbest2(i)-HM2(ceil(rand*HMS),i));
                        Xnew(i)=max(min(Xnew(i),SNPs),1);
                     end 
                else
                     Xnew(i)=HM3(a,i);
                     if rand<PAR                     
                        Xnew(i)=Xnew(i)+((-1)^Xnew(i))*abs(Xbest3(i)-HM3(ceil(rand*HMS),i));
                        Xnew(i)=max(min(Xnew(i),SNPs),1);
                    end 
                end                
          else
                Xnew(i)=ceil(rand*SNPs);
          end
          
          if i==1
              i=i+1;
          elseif i>1 && ~ismember(Xnew(i),Xnew(1:i-1))
              i=i+1;
          end
          
          if ( i-1==dim_epi ) 
              Xtemp=Xnew;
               Xnew=sort(Xnew);
              while ( ismember(Xnew,X,'rows') || ismember(Xnew,X2,'rows') || ismember(Xnew,X3,'rows')|| getMinDistance(Xnew,Center(:,1:dim_epi),dim_epi)<dim_epi-1)
                 j=ceil(rand*dim_epi);
                  r=ceil(rand*SNPs);
                  while ismember(r,Xnew)
                      r=ceil(rand*SNPs);
                  end
                  Xnew(j)=r;
                  Xtemp=Xnew;
                 Xnew=sort(Xnew);
              end
          end          
     end
     
%   Xnew

  [score,score2,score3] = ScoreFunctions(data(:,Xnew),State);

   %%
   Flag = 0;
   
   [fworst,idworst] = max(Fit(:,1));
   [fworst2,idworst2] = max(Fit2(:,2));
   [fworst3,idworst3] = max(Fit3(:,3));
    
        if score<=fworst || score2<=fworst2 || score3<=fworst3
            if  score<=fworst
                Fit(idworst,:)=[score,score2,score3];          
                X(idworst,:)=Xnew; 
                HM(idworst,:)=Xtemp;
                Flag = Flag + 100;
                
            end
            if score2<=fworst2
               Fit2(idworst2,:)=[score,score2,score3];           
               X2(idworst2,:)=Xnew; 
               HM2(idworst2,:)=Xtemp;       
                Flag = Flag + 10;
               
            end
            
            if score3 <= fworst3
               Fit3(idworst3,:)=[score,score2, score3];           
               X3(idworst3,:)=Xnew; 
               HM3(idworst3,:)=Xtemp;       
                Flag = Flag + 1;               
            end

        end
       NC=NC+1; 
 %%  The program is terminted if the Xnew is the solution. 
%  flag = -1;
   if Xnew == CX                   
          canSize1 = canSize1+1;
          Candidate1(canSize1,:) = [CX,score,score2,score3];
          if Flag > 0
              flag = Flag;
          else
              flag = 111;
          end

                       break;
   end
  %% 精英集合管理
  if flag > 0 || mod(NC,maxIterForLocalSearch)==0
                  for i=1:HMS
                      if isempty(Elite1)  %% 同时有2个
                          Elite1 = X(i,:);
                          Efit1 = Fit(i,:);
                          
                          Elite2 = X2(i,:);
                          Efit2 = Fit2(i,:);
                          
                          Elite3 = X3(i,:);
                          Efit3 = Fit3(i,:);
                          
                          
                      else
                              if length(Elite1(:,1))<EliteSize
                                  Elite1=[Elite1;X(i,:)];
                                  Efit1=[Efit1;Fit(i,:)];
                                  
                                  Elite2=[Elite2;X2(i,:)];
                                  Efit2=[Efit2;Fit2(i,:)];
                                  
                                  Elite3=[Elite3;X3(i,:)];
                                  Efit3=[Efit3;Fit3(i,:)];
                              else
                                %%  选择最差的进行替换
                                      [~,eidworst]=max(Efit1(:,1));
                                      [~,eidworst2]=max(Efit2(:,2));
                                      [~,eidworst3]=max(Efit2(:,3));

                                          if Efit1(eidworst,1)> Fit(i,1) %%&& Efit(eidworst,2)> Fit(i,2)
                                              Elite1(eidworst,:)=X(i,:);
                                              Efit1(eidworst,:)=Fit(i,:);
                                          end
                                          if Efit2(eidworst2,2)> Fit(i,2) 
                                              Elite2(eidworst2,:)=X(i,:);
                                              Efit2(eidworst2,:)=Fit(i,:);
                                          end
                                          
                                          if Efit3(eidworst3,3)> Fit(i,3) 
                                              Elite3(eidworst3,:)=X(i,:);
                                              Efit3(eidworst3,:)=Fit(i,:);
                                          end
                                          
                                      
                              end
                      end
                      
                  end
          

   
      %% 获得中心位点
      [~,idebest1]=min(Efit1(:,1));
      [~,idebest2]=min(Efit2(:,2));
      [~,idebest3]=min(Efit3(:,3));
      
      
      E1=Elite1(idebest1,:);
      E2=Elite2(idebest2,:);
    
      E3=Elite3(idebest3,:);
     
          CCC1=Elite1(idebest1,:);
          CCC2=Elite2(idebest2,:);
          CCC3=Elite3(idebest3,:);
          
          
          minDist1=getMinDistance(CCC1,Elite1,dim_epi) ;
           Center=[Center;[CCC1,minDist1]];
           
        if CCC1~=CCC2   
           minDist2=getMinDistance(CCC2,Elite2,dim_epi) ;
           Center=[Center;[CCC2,minDist2]]; 
        end   
        
        
        if CCC3(1,:) ~= CCC1(1,:) & CCC3(1,:) ~= CCC2(1,:)
             minDist3=getMinDistance(CCC3(1,:),Elite3,dim_epi) ;
           Center=[Center;[CCC3(1,:),minDist3]]; 
        end  
        
        
          Alen = min([EliteSize,length(Elite1(:,1)),length(Elite2(:,1)),length(Elite3(:,1))]);
           for i=1:Alen
                      if canSize1==0
                          canSize1=canSize1+1;
                          Candidate1(canSize1,:)=[Elite1(i,:),Efit1(i,:)];
                      else
                          if ~ismember(Elite1(i,:),Candidate1(1:canSize1,1:dim_epi),'rows') 
                              if canSize1<CandidateSize
                                  canSize1=canSize1+1;
                                  Candidate1(canSize1,:)=[Elite1(i,:),Efit1(i,:)];
                              else
                                  [Fitworst,Cind]=max(Candidate1(:,dim_epi+1));  %% dim_epi+1 是标准1， dim_epi+2标准2
                                  if Fitworst>Efit1(i,1) 
                                      Candidate1(Cind,:)=[Elite1(i,:),Efit1(i,:)];
                                  end
                              end
                          end
                      end
                      
                    if canSize2==0
                          canSize2=canSize2+1;
                          Candidate2(canSize2,:)=[Elite2(i,:),Efit2(i,:)];
                    else
                          if ~ismember(Elite2(i,:),Candidate2(1:canSize2,1:dim_epi),'rows') 
                              if canSize2<CandidateSize
                                  canSize2=canSize2+1;
                                  Candidate2(canSize2,:)=[Elite2(i,:),Efit2(i,:)];
                              else
                                  [Fitworst,Cind]=max(Candidate2(:,dim_epi+2));
                                  if Fitworst>Efit2(i,2) 
                                      Candidate2(Cind,:)=[Elite2(i,:),Efit2(i,:)];
                                  end
                              end
                          end
                    end
                      
                     if canSize3==0
                          canSize3=canSize3+1;
                          Candidate3(canSize3,:)=[Elite3(i,:),Efit3(i,:)];
                     else
                        
                          if ~ismember(Elite3(i,:),Candidate3(1:canSize3,1:dim_epi),'rows') 
                              if canSize3<CandidateSize
                                  canSize3=canSize3+1;
                                  Candidate3(canSize3,:)=[Elite3(i,:),Efit3(i,:)];
                              else
                                  [Fitworst,Cind]=max(Candidate3(:,dim_epi+3));
                                  if Fitworst>Efit3(i,3) 
                                      Candidate3(Cind,:)=[Elite3(i,:),Efit3(i,:)];
                                  end
                              end
                          end
                      end
           end
          
%     Candidate(:,1:dim_epi)
    % Center
     
     
    
     
           Elite1=[];
           Efit1=[];
           Elite2=[];
           Efit2=[];
           
           Elite3=[];
           Efit3=[];
           %%重新初始化种群
          %% 小生境排斥
          [X,HM,Fit]=InitHM10(data,HMS,dim_epi,Center);

         X2=X;
         HM2=HM;
         Fit2=Fit;
         
          X3=X;
         HM3=HM;
         Fit3=Fit;
         if ismember(CX,X,'rows')
              canSize1 = canSize1+1;
              Candidate1(canSize1,:) = [CX,score,score2,score3];
              if Flag > 0
                  flag = Flag;
              else
                  flag = 111;
              end

                    
         end
  end

end

Candidate=[Candidate1(1:canSize1,:);Candidate2(1:canSize2,:);Candidate3(1:canSize3,:)];
Candidate = unique(Candidate,'rows');
canSize = length(Candidate(:,1));
totaltime=toc;
% 





