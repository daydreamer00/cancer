function [ pvalues ] = welch_t_test( samp1,samp2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    pvalues = zeros(n1,1);
    for i = 1:n1
        [h, pvalues(i,:)] = ttest2(samp1(i,:),samp2(i,:),'Vartype','unequal');
    end        

end

