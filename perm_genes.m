function [ new_samp1 new_samp2 enum] = perm_genes( samp1,samp2,k )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    enum = combnk(1:n1,k);
    new_n = size(enum,1);
    new_samp1 = single(zeros(new_n,m1));
    new_samp2 = single(zeros(new_n,m2));
   
    w = zeros(new_n,k);
    for i = 1:new_n
        ind = enum(i,:);
        enum_samp1 = samp1(ind,:)';
        enum_samp2 = samp2(ind,:)';
%         w_mat = [w_mat; fish_proj_vec(enum_samp1,enum_samp2)'];
        w(i,:) = fish_proj_vec(enum_samp1,enum_samp2)';
        threshold = - w(i,:)*(mean(enum_samp1)+mean(enum_samp2))'/2;
        
        new_samp1(i,:) = (enum_samp1*w(i,:)')'+repmat(threshold,1,m1);
        new_samp2(i,:) = (enum_samp2*w(i,:)')'+repmat(threshold,1,m2);
    end
    
    mean1 = mean(samp1,2);
    mean2 = mean(samp2,2);
end

