function [ pvalues,enum ] = hotelling_t2_test_perm_genes( samp_pos,samp_neg,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n1 = size(samp_pos,1);
    m1 = size(samp_pos,2);
    n2 = size(samp_neg,1);
    m2 = size(samp_neg,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    
    
    enum = combnk(1:n1,k);
    new_n = size(enum,1);
    pvalues = zeros(new_n,1);
    
    for i = 1:new_n
        if (~mod(i,10000))
            fprintf('iter %d\n',i);
        end
        pos_samp = [samp_pos(enum(i,:),:)];
        neg_samp = [samp_neg(enum(i,:),:)];
%         neg_samp = [samp1(i,npos+1:end); samp2(i,npos+1:end)];
        [t2, pvalues(i,1)] = hotelling_t2_test(pos_samp,neg_samp);
    end

end

