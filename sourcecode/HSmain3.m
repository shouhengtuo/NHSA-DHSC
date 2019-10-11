%% The main program of NHSA-DHSC algorithm

clear;

Dim = 100;
sample_num = 4000;
cacheSize = 20;
epi_dim = 2;
                startIndex = 1; endIndex = 100; 


      pvalue = 0.05/nchoosek(Dim,epi_dim);
      pvalue2 = 1e-4;
       
       HMS = 50;
       CandidateSize =10;
       
      
%         CX = Dim-epi_dim+1:Dim;

       
       folder = 'resultData\';
       dataFile = strcat(folder,'results');
     
        A = {'1st Power','2nd Power', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime','FDR','F1', 'TP', 'TN', 'FP', 'FN', 'Power3'};
    sheet = 1;
   xlRange = 'b1';
   xlswrite(dataFile,A,sheet,xlRange)

  
%% disease model 1 (2-order DME mode)
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.1';
                filepath='modelData\2000CASE_EDM-1_';
                epi_dim = 2;
         
         
%% disease model 2 (3-order DNME mode)
%                 model='Model-2';                
%                 epi_dim = 3;                
%                data = csvread('csvData');  
%                  startIndex = 1;
%                  endIndex = 1;
%%
                
  max_iter = 50000;
   maxIterForLocalSearch = epi_dim*500;              
           % disease loci       
           CX =[Dim - epi_dim + 1 : Dim]
        for dataSetId = startIndex:endIndex
             if dataSetId<10
                noId = strcat('00',num2str(dataSetId));
             elseif dataSetId<100
                noId = strcat('0',num2str(dataSetId));
             else
                noId = num2str(dataSetId);
             end
             data = dlmread(strcat(filepath,noId,'.txt'),'\t',1,0);
        FdataNum = 0; 

        Ac = 0;
        TP = 0;
        FP = 0;
        TN = 0;
        FN = 0;
        
        TP2 = 0;
        FP2 = 0;
        TN2 = 0;
        FN2 = 0;
         power1 = 0;  
         power2 = 0; 
         power3 = 0;
         power4 = 0;
         K2_power = 0;
         Gini_power = 0;
         JE_power = 0;
        Evalutation_Times = [];
        TIME = [];
         succ = 0; 
        succEvalutation_Times = [];
        succTIME = [];
        
       
         
          
%% Search k-snp loci using Harmony search algorithm 
                [Candidate,canSize,Nc,runtime,flag] = NHSA3(data,epi_dim,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX);
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
                 fprintf(2,'\n success search time(%f)|(%d),   success ',  runtime, Nc);
            else 
                fprintf('\n   search time(%f)|(%d)  **** fail! *** ',runtime,Nc );
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
    
    
   Results=[Power1,Power2, K2_Power, Gini_Power,JE_Power, TPR,SPC, PPV, ACC, E_Times, RunTime, succE_Times, succRunTime, FDR,F1, TP, TN, FP, FN, Power3];
    
    
   % A = {'Power1','Power2','Power3', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime'};
    sheet = 1;
   xlRange = 'B2' ;
   xlswrite(dataFile,Results,sheet,xlRange)




