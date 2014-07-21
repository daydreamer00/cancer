function [t2_statistic,p_value] = hotelling_t2_test(Sample_1,Sample_2)
%========================================================
% function t=hotelling_t2_test_t_test(Sample_1,Sample_2)
%--------------------------------------------------------
% [Instruction]: two-samples' hotelling t2 test
% Version 1.0, by Percy 
% Date: 2014-5-6
%========================================================
if size(Sample_1,1)~=size(Sample_2,1)
    error('Wrong input!');
else
    m = size(Sample_1,1);
end
n1_rev = size(Sample_1,2);
n2_rev = size(Sample_2,2);
M1 = mean(Sample_1,2);
M2 = mean(Sample_2,2);
SS1 = cov(Sample_1');
SS2 = cov(Sample_2');
SS = (SS1*(n1_rev-1)+SS2*(n2_rev-1))./(n1_rev+n2_rev-2);
t2_statistic = (M1-M2)'*pinv(SS)*(M1-M2)*n1_rev*n2_rev/(n1_rev+n2_rev);
t2_statistic = (n1_rev+n2_rev-m-1)/(m*(n1_rev+n2_rev-2))*t2_statistic;
df1 = m;
df2 = n1_rev+n2_rev-m-1;
% p_value = 1-cdf('f',t2_statistic,df1,df2);   % p-value
p_value = cdf('f',1/t2_statistic,df2,df1);
end