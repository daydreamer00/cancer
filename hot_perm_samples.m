function [ pvalues ] = hot_perm_samples( samp1,samp2,iter_n,n_pos )
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
    
    [ ori_hot_pvalues ] = hotelling_t2_test_batch( samp1,samp2,n_pos );
    
    cnt = zeros(n1,1);
    
    for i = 1:iter_n
        fprintf('iter %d start\n',i);
        index = randperm(m1);
        pos_i = index<=n_pos;
        neg_i = index>n_pos;
        
        samp1_pos = samp1(:,pos_i);
        samp1_neg = samp1(:,neg_i);
        samp2_pos = samp2(:,pos_i);
        samp2_neg = samp2(:,neg_i);        
        
        new_samp1 = [samp1_pos samp1_neg];
        new_samp2 = [samp2_pos samp2_neg];
        
        [ perm_hot_pvalues ] = hotelling_t2_test_batch( new_samp1,new_samp2,n_pos );
        
        cnt = cnt+(perm_hot_pvalues<ori_hot_pvalues);      
        
    end
    N = iter_n;

    pvalues = cnt/N;

end

