function [ res_table ] = compute_acc_perm_genes( samp_pos,samp_neg,k,iter)
%UNTITLED2 Summary of this function goes here
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
    res_table = zeros(new_n,3);
    
     train_r = 0.7;
     n_pos = m1;
     n_neg = m2;
     
     n_train_samp_pos = fix(n_pos*train_r);
     n_train_samp_neg = fix(n_neg*train_r);
     n_test_samp_pos = n_pos - n_train_samp_pos;
     n_test_samp_neg = n_neg - n_train_samp_neg;
     
     
     
     labels = [ones(m1,1);zeros(m2,1)];
     
     for i = 1:new_n
         if (~mod(i,10000))
             fprintf('iter %d\n',i);
         end
         pos_samp = [samp_pos(enum(i,:),:)];
         neg_samp = [samp_neg(enum(i,:),:)];
         %         neg_samp = [samp1(i,npos+1:end); samp2(i,npos+1:end)];
         
         cum_acc_samp_pos_train = 0;
         cum_acc_samp_neg_train = 0;
         cum_acc_samp_pos_test = 0;
         cum_acc_samp_neg_test = 0;
         
         train_eval = zeros(1,7);
         test_eval = zeros(1,7);
         
         for j = 1:iter
             index_samp_pos = randperm(n_pos);
             index_samp_neg = randperm(n_neg);
             
             train_i_pos = index_samp_pos(1:n_train_samp_pos);
             train_i_neg = index_samp_neg(1:n_train_samp_neg);
             test_i_pos = index_samp_pos(n_train_samp_pos+1:end);
             test_i_neg = index_samp_neg(n_train_samp_neg+1:end);
             
             samp_pos_train = pos_samp(:,train_i_pos);
             samp_neg_train = neg_samp(:,train_i_neg);
             samp_train = [samp_pos_train samp_neg_train];
             
             samp_pos_test = pos_samp(:,test_i_pos);
             samp_neg_test = neg_samp(:,test_i_neg);
             samp_test = [samp_pos_test samp_neg_test];
             
             w = fish_proj_vec(samp_pos_train',samp_neg_train')';
             
             thres = - w*(mean(samp_pos_train,2)+mean(samp_neg_train,2))/2;
             
             predict_train_pos_label = (w*samp_pos_train+thres)>0;
             predict_train_neg_label = (w*samp_neg_train+thres)>0;
             predict_test_pos_label = (w*samp_pos_test+thres)>0;
             predict_test_neg_label = (w*samp_neg_test+thres)>0;
             
             true_train_label = [ones(n_train_samp_pos,1);zeros(n_train_samp_neg,1)];
             true_test_label = [ones(n_test_samp_pos,1); zeros(n_test_samp_neg,1)];
             
             train_labels = [predict_train_pos_label predict_train_neg_label];
             test_labels = [predict_test_pos_label predict_test_neg_label];
             
             train_eval = train_eval + Evaluate(true_train_label,train_labels');
             test_eval = test_eval + Evaluate(true_test_label,test_labels');
                         
                          
             acc_samp_pos_train = sum(predict_train_pos_label==1)/n_train_samp_pos;
             acc_samp_neg_train = sum(predict_train_neg_label==0)/n_train_samp_neg;
             acc_samp_pos_test = sum(predict_test_pos_label==1)/n_test_samp_pos;
             acc_samp_neg_test = sum(predict_test_neg_label==0)/n_test_samp_neg;
             
             cum_acc_samp_pos_train = cum_acc_samp_pos_train+acc_samp_pos_train;
             cum_acc_samp_neg_train = cum_acc_samp_neg_train+acc_samp_neg_train;
             cum_acc_samp_pos_test = cum_acc_samp_pos_test+acc_samp_pos_test;
             cum_acc_samp_neg_test = cum_acc_samp_neg_test+acc_samp_neg_test;
         end
         train_eval = train_eval/iter;
         test_eval = test_eval/iter;
         acc_train = (cum_acc_samp_pos_train+cum_acc_samp_neg_train)/2/iter;
         acc_test = (cum_acc_samp_pos_test+cum_acc_samp_neg_test)/2/iter;
         res_table(i,:) = [acc_train acc_test (acc_train+acc_test)/2];
     end
     res_table = [enum res_table];     
end

