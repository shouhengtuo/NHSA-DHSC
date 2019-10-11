function [P_value, P_value2] = permutation(snp_com,state,permutation_times,P_value0, P_value02)
%     P_value0 = Gtest_score(snp_com,state);
    pvalue = zeros(permutation_times,1);
    L = length(state);
    for k = 1:permutation_times
        A = randperm(L);
        Astate = state(A);
        %state = round(rand(L,1));
        pvalue(k) = Gtest_score(snp_com,Astate);
    end
    
    P_value = sum(pvalue <= P_value0)/permutation_times;
    P_value2 = sum(pvalue <= P_value02)/permutation_times;
end