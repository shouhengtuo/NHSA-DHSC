
 clear;

Dim = 100;
sample_num = 3000;
cacheSize = 20;
pvalue2 = 0.01;

%% harmony search algorithm parameters setting
       HMS = 100;
       CandidateSize = 100;

%%
folder = 'NHSA_NDMEData2019\';

dataFile = strcat(folder,'MultiPower',num2str(Dim));
    A = {'1st Power','2nd Power','3rd Power', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime','FDR','F1'};
    sheet = 1;
   xlRange = 'b1';
   xlswrite(dataFile,A,sheet,xlRange)
   
for dataIndex = [1]
        switch(dataIndex)
            case 1
                epi_dim = 3;
                startIndex = 1; endIndex = 100;
                filepath = '.\threewayBests\'; filename = 'threewayBests'; 
            case 2
                epi_dim = 3;
                startIndex = 401; endIndex = 500;
                filepath = '.\HWthreewayBests\';filename = 'HWthreewayBests'; 
            case 3
                epi_dim = 4;
                startIndex = 101; endIndex = 200; 
                filepath = '.\fourwayBests\';filename = 'fourwayBests';
            case 4
                epi_dim = 4;
                startIndex = 1101; endIndex = 1200;
                filepath = '.\fourwayNoLowBests\';    filename = 'fourwayNoLowBests';
            case 5
                epi_dim = 4;
                startIndex = 501; endIndex = 600;
                filepath = '.\HWfourwayBests\';  filename = 'HWfourwayBests';
            case 6
                epi_dim = 5;
                startIndex = 201; endIndex = 300;
                filepath = '.\fivewayBests\';  filename = 'fivewayBests';
            case 7 
               epi_dim = 5;
                startIndex = 1201; endIndex = 1300;
                filepath = '.\fivewayNoLowBests\';  filename = 'fivewayNoLowBests';
           case 8 
               epi_dim = 5;
                startIndex = 601; endIndex = 700;
                filepath = '.\HWfivewayBests\';       filename = 'HWfivewayBests'; 
        end
       %% harmony search algorithm parameters setting  
       max_iter = 3000 * epi_dim^3;
       maxIterForLocalSearch = max_iter/5;
      
       HMS = 100;
       CandidateSize = 10;
       
        CX = Dim-epi_dim+1:Dim;

        pvalue = 0.05/nchoosek(Dim,epi_dim);

        Ac = 0;
         TP = 0;
        FP = 0;
        TN = 0;
        FN = 0;
         power1 = 0;  
         power2 = 0; 
         power3 = 0;
         K2_power = 0;
         Gini_power = 0;
         JE_power = 0;
        Evalutation_Times = [];
        TIME = [];
         succ = 0; 
        succEvalutation_Times = [];
        succTIME = [];
        
        for dataSet = startIndex:endIndex
           % a = dlmread(strcat(filepath,strcat('best',num2str(dataSet),'.txt')),',',1,0);
             a = dlmread(strcat(filepath,strcat('best',num2str(dataSet),'.txt')),'\t',1,0);
            AA = 0; Aa = 0; aa = 0;
                for i =1:1500%1:sample_num
                    for j = 1:epi_dim
                        if a(i,j) == 2
                            aa = aa + 1;
                        elseif a(i,j) == 1
                            Aa = Aa + 1;
                        else
                            AA = AA + 1;
                        end
                    end
                end
                AA = 2*AA /(sample_num*epi_dim);
                Aa = 2*Aa /(sample_num*epi_dim);
                aa = 2*aa /(sample_num*epi_dim);
                
                b = zeros(sample_num,Dim-epi_dim);
                for i = 1 : sample_num
                    for j = 1:Dim-epi_dim
                        r = rand;
                        if r <= aa
                            b(i,j) = 2;
                        elseif r <= aa+Aa
                            b(i,j) = 1;
                        else
                            b(i,j) = 0;
                        end
                          % b(i,j) = fix(rand*3);
                    end
                end

            data = [b,a];
            
    %% Search k-snp loci using Harmony search algorithm 
                [Candidate,canSize,Nc,runtime,flag] = HS_FOR_multiLOCI12(data,epi_dim,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX);
                     Evalutation_Times = [Evalutation_Times , Nc];
             TIME = [TIME  runtime];
             if flag == 111
                 K2_power = K2_power + 1;
                 Gini_power = Gini_power + 1;
                JE_power = JE_power + 1;
             elseif flag == 110
                  K2_power = K2_power + 1;
                 Gini_power = Gini_power + 1;
             elseif flag == 101
                 K2_power = K2_power + 1;                
                JE_power = JE_power + 1;
            elseif flag == 11
                JE_power = JE_power + 1;
                 Gini_power = Gini_power + 1;
            elseif flag == 10
                Gini_power = Gini_power + 1;
            elseif flag == 1
                JE_power = JE_power + 1;
             end
                 
                 
         % end harmony search  ****************************************
              [SNP_COM1, bestId1] = sort(Candidate(:,epi_dim+1));
              [SNP_COM2, bestId2] = sort(Candidate(:,epi_dim+2));
              [SNP_COM3, bestId3] = sort(Candidate(:,epi_dim+3));
              
           
           %% 1st stage power  
            if flag > 0
                power1 = power1 + 1;
                 succ = succ + 1;
                succEvalutation_Times = [succEvalutation_Times, Nc];
                 succTIME = [succTIME  runtime];
                 fprintf(2,'\ndataSet %3d:success search time(%f)|(%d),   success %d/%d ',dataSet, runtime, Nc,succ,dataSet-startIndex+1);
            else 
                fprintf('\ndataSet %3d: search time(%f)|(%d), %d/%d**** fail! *** ',dataSet,runtime,Nc,succ,dataSet-startIndex+1);
            end

        %%2nd stage: G-test
       G_set = [];
            for i = 1:canSize
                P_value = Gtest_score(Candidate(i,1:epi_dim),data(:,Dim+1));
               if P_value < pvalue
                    G_set = [G_set;Candidate(i,1:epi_dim)];
               end
            end
          
           if ~isempty(G_set) && ismember(CX,G_set,'rows')
              power2 = power2 + 1;
              G_setSize = length(G_set(:,1));
              fprintf('G-set size = %d\n',  G_setSize);
           else
                G_setSize = 0;
               fprintf('G-set size = %d\n',  G_setSize);
           end
%% 3rd stage 
           %% 3rd stage 
                tp = 0;
                fp = 0;
                fn = 0;
                tn = 0;
            for i = 1:  G_setSize              
                   perm_Pvalue = permutation(data(:,G_set(i,:)),data(:,Dim+1),100,pvalue,pvalue2);

                  if isequal(G_set(i,:), Dim-epi_dim+1:Dim )
                       if perm_Pvalue < pvalue
                             tp = tp + 1;  
                       else
                            fn = fn + 1;
                      end
                  else
                        if perm_Pvalue < pvalue
                           fp = fp + 1;
                       else
                           tn = tn + 1;
                        end
                  end
                      
            end
                     
                    if tp > 0 
                        TP = TP + 1;
                        power3 = power3 + 1;
                    else
                        FN = FN + 1;
                    end
                        

            

                      if tn + fp == 0
                          TN = TN + 1;
                      elseif fp > 0
                          FP = FP + 1;
                      else
                          TN = TN + 1;
                      end
         
              
              
        end
    
        
  Power1 = power1/100;
   Power2 = power2/100;
   Power3 = power3/100;
   K2_Power = K2_power/100;
   Gini_Power = Gini_power/100;
   JE_Power = JE_power/100;
   E_Times = median( Evalutation_Times);
   RunTime = median( TIME);
   succE_Times = median (succEvalutation_Times);
   succRunTime = median( succTIME);
   
   
  TPR = TP/(TP + FN);
  SPC = TN/(FP + TN);
  PPV = TP/(TP + FP);
  ACC = (TP + TN)/(TP + TN + FN + FP);
  FDR = 1 - PPV;
  F1 = 2*TP/(2*TP + FP + FN);

   
    dataFile = strcat(folder,'MultiPower', num2str(Dim));
    Results=[Power1,Power2,Power3, K2_Power, Gini_Power,JE_Power, TPR,SPC, PPV, ACC, E_Times, RunTime, succE_Times, succRunTime, FDR,F1];
    
    
   % A = {'Power1','Power2','Power3', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime'};
    sheet = 1;
   xlRange = strcat('B', num2str(dataIndex+1), ': Q', num2str(dataIndex+1)) ;
   xlswrite(dataFile,Results,sheet,xlRange)
    sheet = 2; % evaluation times
   xlRange = strcat('B', num2str(dataIndex+1));
   xlswrite(dataFile,Evalutation_Times,sheet,xlRange) 
   
   sheet = 3; % evaluation times
   xlRange = strcat('B', num2str(dataIndex+1));
   xlswrite(dataFile,TIME,sheet,xlRange)    
end



