%% The main program of NHSA-DHSC algorithm

clear;

Dim = 100;
sample_num = 4000;
cacheSize = 20;
epi_dim = 2;
                startIndex = 1; endIndex = 100; 


      pvalue = 0.05/nchoosek(Dim,epi_dim);
      pvalue2 = 1e-4;
       max_iter = 50000;
       maxIterForLocalSearch = epi_dim*500;
       HMS = 50;
       CandidateSize =10;
       
      
%         CX = Dim-epi_dim+1:Dim;

       
       folder = 'resultData\';
       dataFile = strcat(folder,'results');
     
        A = {'1st Power','2nd Power', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime','FDR','F1', 'TP', 'TN', 'FP', 'FN', 'Power3'};
    sheet = 1;
   xlRange = 'b1';
   xlswrite(dataFile,A,sheet,xlRange)

  
%% data model 
                model='Model-2';parameter='H2=0.02,PD=0.1,MAF=0.1';
                filepath='modelData\2000CASE_EDM-1_';
          
 %% disease loci       
           CX =[Dim-1, Dim];

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
        
         for dataSetId = startIndex:endIndex
             if dataSetId<10
                noId = strcat('00',num2str(dataSetId));
             elseif dataSetId<100
                noId = strcat('0',num2str(dataSetId));
             else
                noId = num2str(dataSetId);
             end
             data = dlmread(strcat(filepath,noId,'.txt'),'\t',1,0);
             
             
           CX =  [Dim - epi_dim + 1 : Dim];
  %% swap the diease-causing loci    randomly, which aims to test the ability for detecting different SNP-combination, many exist search algorithms have the preference to the position of diease loci.
             rr = ceil(rand(1,2)*Dim);
                   while rr(1) == rr(2)
                        rr = ceil(rand(1,2)*Dim);
                   end
             temp = data(:,rr);
             data(:,rr) = data(:,CX );
             data(:,CX) = temp;
             
             CX = sort( rr)           
                   
             
             
             
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
             elseif flag == 100 
                 K2_power = K2_power + 1;
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
                 fprintf('\ndataSet %3d:success search time(%f),   success    ',dataSetId, runtime);
            else 
                fprintf('\ndataSet %3d: search time(%f) **** fail! *** ',dataSetId,runtime);
            end

            %% 2nd stage: G-test 
                   G_set = [];
                        for i = 1:canSize
                            P_value = Gtest_score(data(:,Candidate(i,1:epi_dim)),data(:,Dim+1));
                           if P_value < pvalue
                                G_set = [G_set;Candidate(i,1:epi_dim)];
                           end
                        end

                          TN = TN + canSize;

                       if isempty(G_set)
                           FN = FN + 1;

                       elseif ismember(CX,G_set,'rows')

                          power2 = power2 + 1;
                          TP = TP + 1;
                          G_setSize =  length(G_set(:,1));
                          TN = TN - G_setSize;
                       else
                            G_setSize =  length(G_set(:,1));
                            FP = FP + G_setSize;
                           TN = TN - G_setSize;
                       end

%         %%2nd stage: chi-square test
%                    G_set = [];
%                         for i = 1:canSize
%                             [P_value, Gvalue] = Chi_square2(data(:,Candidate(i,1:epi_dim)),data(:,Dim+1));
%                            if P_value < pvalue
%                                 G_set = [G_set;Candidate(i,1:epi_dim)];
%                            end
%                         end
% 
%                           TN2 = TN2 + canSize;
% 
%                        if isempty(G_set)
%                            FN2 = FN2+ 1;
% 
%                        elseif ismember(CX,G_set,'rows')
% 
%                           power3 = power3 + 1;
%                           TP2 = TP2 + 1;
%                           G_setSize =  length(G_set(:,1));
%                           TN2 = TN2 - G_setSize;
%                        else
%                             G_setSize =  length(G_set(:,1));
%                             FP2 = FP2 + G_setSize;
%                            TN2 = TN2- G_setSize;
%                        end    
%               
              
         end
    
        


    
  Power1 = power1/100;
   Power2 = power2/100;

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
  
  Power3 = power3/100;
  
%   TPR2 = TP2/(TP2 + FN2);
%   SPC2 = TN2/(FP2 + TN2);
%   PPV2 = TP2/(TP2 + FP2);
%   ACC2 = (TP2 + TN2)/(TP2 + TN2 + FN2 + FP2);
%   FDR2 = 1 - PPV2;
%   F12 = 2*TP2/(2*TP2 + FP2 + FN2);

   
   
    Results=[Power1,Power2, K2_Power, Gini_Power,JE_Power, TPR,SPC, PPV, ACC, E_Times, RunTime, succE_Times, succRunTime, FDR,F1, TP, TN, FP, FN, Power3];
    
    
   % A = {'Power1','Power2','Power3', 'K2_Power', 'Gini_Power','JE_Power', 'TPR','SPC', 'PPV', 'ACC', 'E_Times', 'RunTime', 'succE_Times', 'succRunTime'};
    sheet = 1;
   xlRange = 'B2' ;
   xlswrite(dataFile,Results,sheet,xlRange)
    
     
