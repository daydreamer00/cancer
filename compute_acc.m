function [ acc_train, acc_test ] = compute_acc( samp1, samp2, iter )
%COMPUTE_ACC compute the accuraccy of each gene(combinantion)
%   Detailed explanation goes here
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    
    train_r = 0.7;
    
    n_train_samp1 = fix(m1*train_r);
    n_train_samp2 = fix(m2*train_r);
    n_test_samp1 = m1 - n_train_samp1;
    n_test_samp2 = m2 - n_train_samp2;
    
    cum_acc_samp1_train = zeros(n1,1);
    cum_acc_samp2_train = zeros(n1,1);
    cum_acc_samp1_test = zeros(n1,1);
    cum_acc_samp2_test = zeros(n1,1);
    
    for i = 1:iter
%         fprintf('iter %d start\n',i);
        index_samp1 = randperm(m1);
        index_samp2 = randperm(m2);
        samp1_train = samp1(:,index_samp1(1:n_train_samp1));
        samp2_train = samp2(:,index_samp2(1:n_train_samp2));

        samp1_mean_train = mean(samp1_train,2);
        samp2_mean_train = mean(samp2_train,2);
        
        threshold = (samp1_mean_train+samp2_mean_train)/2;
        
        acc_samp1_train = sum(samp1_train>=repmat(threshold,1,n_train_samp1),2)/n_train_samp1;
        acc_samp2_train = sum(samp2_train<repmat(threshold,1,n_train_samp2),2)/n_train_samp2;
        acc_samp1_test = sum(samp1(:,index_samp1(n_train_samp1+1:end))>=repmat(threshold,1,n_test_samp1),2)/n_test_samp1;
        acc_samp2_test = sum(samp2(:,index_samp2(n_train_samp2+1:end))<repmat(threshold,1,n_test_samp2),2)/n_test_samp2;
        
        cum_acc_samp1_train = cum_acc_samp1_train+acc_samp1_train;
        cum_acc_samp2_train = cum_acc_samp2_train+acc_samp2_train;
        cum_acc_samp1_test = cum_acc_samp1_test+acc_samp1_test;
        cum_acc_samp2_test = cum_acc_samp2_test+acc_samp2_test;
    end     
    acc_train = (cum_acc_samp1_train+cum_acc_samp2_train)/2/iter;
    acc_test = (cum_acc_samp1_test+cum_acc_samp2_test)/2/iter;
end

