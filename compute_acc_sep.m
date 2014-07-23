function [ acc_train, acc_test,train_eval, test_eval ] = compute_acc_sep( samp1,samp2,k_neigh,projtype,disttype,meantype,n_pos,iter )
%COMPUTE_ACC compute the accuraccy of each gene(combinantion)
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
    n_neg = m1 - n_pos;
    
    train_r = 0.7;
    
    n_train_samp_pos = fix(n_pos*train_r);
    n_train_samp_neg = fix(n_neg*train_r);
    n_test_samp_pos = n_pos - n_train_samp_pos;
    n_test_samp_neg = n_neg - n_train_samp_neg;
    
    cum_acc_samp_pos_train = zeros(n1,1);
    cum_acc_samp_neg_train = zeros(n1,1);
    cum_acc_samp_pos_test = zeros(n1,1);
    cum_acc_samp_neg_test = zeros(n1,1);
    
    train_eval = zeros(n1,7);
    test_eval = zeros(n1,7);
    
    labels = [zeros(m1,1);ones(m2,1)];
    
    for i = 1:iter
        fprintf('iter %d start\n',i);
        index_samp_pos = randperm(n_pos);
        index_samp_neg = randperm(n_neg);
        
        train_i_pos = index_samp_pos(1:n_train_samp_pos);
        train_i_neg = index_samp_neg(1:n_train_samp_neg);
        test_i_pos = index_samp_pos(n_train_samp_pos+1:end);
        test_i_neg = index_samp_neg(n_train_samp_neg+1:end);
        
        samp1_pos_train = samp1(:,train_i_pos);
        samp1_neg_train = samp1(:,n_pos+train_i_neg);
        samp1_train = [samp1_pos_train samp1_neg_train];
        
        samp2_pos_train = samp2(:,train_i_pos);
        samp2_neg_train = samp2(:,n_pos+train_i_neg);
        samp2_train = [samp2_pos_train samp2_neg_train];
        
        samp1_pos_test = samp1(:,test_i_pos);
        samp1_neg_test = samp1(:,n_pos+test_i_neg);
        samp1_test = [samp1_pos_test samp1_neg_test];
        
        samp2_pos_test = samp2(:,test_i_pos);
        samp2_neg_test = samp2(:,n_pos+test_i_neg);
        samp2_test = [samp2_pos_test samp2_neg_test];
        
        [ new_samp_pos_train, new_samp_neg_train ,w,thres] = k_near_proj( samp1_train,samp2_train,k_neigh,projtype,disttype,meantype,n_train_samp_pos,0,[]);
        
        new_samp_pos_test = zeros(n1,n_test_samp_pos);
        new_samp_neg_test = zeros(n1,n_test_samp_neg);
        
%         mean_samp1_test = mean(samp1_test,2);
%         mean_samp2_test = mean(samp2_test,2);
        
        mean_samp1_pos_test = mean(samp1_pos_test,2);
        mean_samp1_neg_test = mean(samp1_neg_test,2);
        mean_samp2_pos_test = mean(samp2_pos_test,2);
        mean_samp2_neg_test = mean(samp2_neg_test,2);
        
%         mean_samp1_test = (mean_samp1_pos_test*n_test_samp_pos+mean_samp1_neg_test*n_test_samp_neg)/(n_test_samp_neg+n_test_samp_pos);
%         mean_samp2_test = (mean_samp2_pos_test*n_test_samp_pos+mean_samp2_neg_test*n_test_samp_neg)/(n_test_samp_neg+n_test_samp_pos);
       
        mean_samp1_test = (mean_samp1_pos_test+mean_samp1_neg_test)/2;
        mean_samp2_test = (mean_samp2_pos_test+mean_samp2_neg_test)/2;

        for j = 1:n1
%             thres(j,:) = -w(j,:)*[mean_samp1_test(j); mean_samp2_test(j);];
            new_samp_pos_test(j,:) = w(j,:)*[samp1_pos_test(j,:); samp2_pos_test(j,:)]+repmat(thres(j,:),1,n_test_samp_pos);
            new_samp_neg_test(j,:) = w(j,:)*[samp1_neg_test(j,:); samp2_neg_test(j,:)]+repmat(thres(j,:),1,n_test_samp_neg);
        end
        
%         for j = 1:n1
%             new_samp_pos_test(j,:) = w(j,:)*[samp1_pos_test(j,:); samp2_pos_test(j,:)]+repmat(thres(j,:),1,n_test_samp_pos);
%             new_samp_neg_test(j,:) = w(j,:)*[samp1_neg_test(j,:); samp2_neg_test(j,:)]+repmat(thres(j,:),1,n_test_samp_neg);
%         end
            
        train_pos_labels = new_samp_pos_train<zeros(n1,n_train_samp_pos);
        train_neg_labels = new_samp_neg_train<zeros(n1,n_train_samp_neg);
        test_pos_labels = new_samp_pos_test<zeros(n1,n_test_samp_pos);
        test_neg_labels = new_samp_neg_test<zeros(n1,n_test_samp_neg);
        
        true_train_label = [zeros(n_train_samp_pos,1);ones(n_train_samp_neg,1)];
        true_test_label = [zeros(n_test_samp_pos,1); ones(n_test_samp_neg,1)];
        
        train_labels = [train_pos_labels train_neg_labels];
        test_labels = [test_pos_labels test_neg_labels];
        
        for j = 1:n1
            train_eval(j,:) = train_eval(j,:) + Evaluate(true_train_label,train_labels(j,:)');
            test_eval(j,:) = test_eval(j,:) + Evaluate(true_test_label,test_labels(j,:)');
        end
        
        acc_samp_pos_train = sum(new_samp_pos_train>=zeros(n1,n_train_samp_pos),2)/n_train_samp_pos;
        acc_samp_neg_train = sum(new_samp_neg_train<zeros(n1,n_train_samp_neg),2)/n_train_samp_neg;
        acc_samp_pos_test = sum(new_samp_pos_test>=zeros(n1,n_test_samp_pos),2)/n_test_samp_pos;
        acc_samp_neg_test = sum(new_samp_neg_test<zeros(n1,n_test_samp_neg),2)/n_test_samp_neg;
        
        cum_acc_samp_pos_train = cum_acc_samp_pos_train+acc_samp_pos_train;
        cum_acc_samp_neg_train = cum_acc_samp_neg_train+acc_samp_neg_train;
        cum_acc_samp_pos_test = cum_acc_samp_pos_test+acc_samp_pos_test;
        cum_acc_samp_neg_test = cum_acc_samp_neg_test+acc_samp_neg_test;
    end
    train_eval = train_eval/iter;
    test_eval = test_eval/iter;
    acc_train = (cum_acc_samp_pos_train+cum_acc_samp_neg_train)/2/iter;
    acc_test = (cum_acc_samp_pos_test+cum_acc_samp_neg_test)/2/iter;
end

