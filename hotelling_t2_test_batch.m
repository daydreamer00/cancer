function [ pvalues ] = hotelling_t2_test_batch( samp1,samp2,npos)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    if(m1~=m2)
        error('not same number of samples');
    end
    pvalues = zeros(n1,1);
    
    for i = 1:n1
        pos_samp = [samp1(i,1:npos); samp2(i,1:npos)];
        neg_samp = [samp1(i,npos+1:end); samp2(i,npos+1:end)];
        [t2, pvalues(i,1)] = hotelling_t2_test(pos_samp,neg_samp);
    end

end

