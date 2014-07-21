%load data
ori_data = xlsread('miRNA_with_metastasis.xlsx');
save('ori_data.mat','ori_data');
load('ori_data.mat');
samp1 = ori_data(5:end,ori_data(1,:)==1);
samp2 = ori_data(5:end,ori_data(1,:)==2);

n1 = size(samp1,1);

n_pos = sum(ori_data(4,:)==0)/2;

%generate random data
% samp1 = randn(size(samp1));
% samp2 = randn(size(samp1));

k_neigh = 1
proj_type = 2 % 1:fisher 2:pSVM 3:libSVM 4: tuned libSVM
dist_type = 1
mean_type = 1 % 1:normal 2:weighted

% t = 180;
% samp1(1,:) = samp1(t,:);
% samp2(1,:) = samp2(t,:);

[ new_samp_pos, new_samp_neg ] = k_near_proj( samp1,samp2,k_neigh,proj_type,dist_type,mean_type,n_pos,1,[]);
new_samp = [new_samp_pos new_samp_neg];
if(proj_type ==1 )
    save(sprintf('%d nn fisher.txt',k_neigh),'new_samp','-ascii');
elseif(proj_type==2)
    save(sprintf('%d nn svm.txt',k_neigh),'new_samp','-ascii');
elseif(proj_type==3)
    save(sprintf('%d nn libsvm.txt',k_neigh),'new_samp','-ascii');
end

[ ttest_pvalues ] = welch_t_test( new_samp_pos,new_samp_neg );
[minpvalue min_i] = min(ttest_pvalues);
min_i
figure(1);
plotHist2(new_samp_pos(min_i,:),new_samp_neg(min_i,:));
[sorted_p, sort_i] = sort(ttest_pvalues);
sorted_p_ttest = [sorted_p sort_i];

% figure(1);
% i=216
% plotHist2(new_samp_pos(i,:),new_samp_neg(i,:));

k = 1;
[ new_samp1, new_samp2, enum] = perm_genes(  new_samp_pos,new_samp_neg,k );

% %permutation t test for projected
% [perm_pvalues ] = perm_samples( new_samp1,new_samp2,5000,1); 
% [sorted_p, sort_i] = sort(perm_pvalues);
% sorted_p_perm = [sorted_p sort_i];

%permutation hotelling
[perm_pvalues ] = hot_perm_samples( new_samp1,new_samp2,5000,1); 
[sorted_p, sort_i] = sort(perm_pvalues);
sorted_p_perm = [sorted_p sort_i];

[ acc_train, acc_test ,train_eval, test_eval] = compute_acc_sep( samp1,samp2,k_neigh,proj_type,dist_type,mean_type,n_pos,500 );
avg_acc = (acc_train+acc_test)/2;
acc_sep = [acc_train acc_test avg_acc];
[ acc_train, acc_test ] = compute_acc( new_samp1,new_samp2,5000 );
avg_acc = (acc_train+acc_test)/2;
acc = [acc_train acc_test avg_acc];

%hotelling
[ hot_pvalues ] = hotelling_t2_test_batch( samp1,samp2,n_pos );

res_table_perm = [[1:n1]' hot_perm_pvalues hot_pvalues ttest_pvalues acc_sep train_eval test_eval acc];

% samp1_pos = samp1(:,1:n_pos);
% samp1_neg = samp1(:,n_pos+1:end);
% samp2_pos = samp2(:,1:n_pos) - samp1_pos;
% samp2_neg = samp2(:,n_pos+1:end)-samp1_neg;
% 
% [ hot_pvalues,enum ] = hotelling_t2_test_perm_genes( samp2_pos,samp2_neg,3);
% hot_pvalues = [enum hot_pvalues];
% save('perm_hot_pvalue_cancerMinus.txt','hot_pvalues','-ascii');


% 
% %only use cancer
% samp2_pos = samp2(:,1:n_pos);
% samp2_neg = samp2(:,n_pos+1:end);
% [ pvalues ] = welch_t_test( samp2_pos,samp2_neg );
% [minpvalue min_i] = min(pvalues);
% min_i
% figure(2);
% plotHist2(samp2_pos(min_i,:),samp2_neg(min_i,:));
% [sorted_p, sort_i] = sort(pvalues);
% sorted_p_cancer = [sorted_p sort_i];
% 
% %only use normal
% samp1_pos = samp1(:,1:n_pos);
% samp1_neg = samp1(:,n_pos+1:end);
% [ pvalues ] = welch_t_test( samp1_pos,samp1_neg );
% [minpvalue min_i] = min(pvalues);
% min_i
% figure(3);
% plotHist2(samp1_pos(min_i,:),samp1_neg(min_i,:));
% [sorted_p, sort_i] = sort(pvalues);
% sorted_p_norm = [sorted_p sort_i];
% 
% %plot gene i
% i=279
% %projected
% figure(1);
% plotHist2(new_samp_pos(i,:),new_samp_neg(i,:));
% %original
% figure(2);
% plotHist2(samp2_pos(i,:),samp2_neg(i,:));
% figure(3);
% plotHist2(samp1_pos(i,:),samp1_neg(i,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %train SVM with all genes
% ori_data = xlsread('miRNA_with_metastasis.xlsx');
% samp1 = ori_data(5:end,ori_data(1,:)==1);
% samp2 = ori_data(5:end,ori_data(1,:)==2);
% n_pos = sum(ori_data(4,:)==0)/2;
% n_neg = size(samp1,2) - n_pos;
% train_r = 0.7;
% n_train_pos = fix(n_pos*train_r);
% n_train_neg = fix(n_neg*train_r);
% train_label = [ ones(n_train_pos,1); zeros(n_train_neg,1)];
% n_test_pos = n_pos - n_train_pos;
% n_test_neg = n_neg - n_train_neg;
% test_label = [ ones(n_test_pos,1); zeros(n_test_neg,1)];
% samp = samp1';
% train_mat = samp([1:n_train_pos n_pos+(1:n_train_neg)],:);
% test_mat = samp([n_train_pos+(1:n_test_pos) n_pos+n_train_neg+(1:n_test_neg)],:);
% [stand_train_mat mean_a std_a ]  = standardize(train_mat);
% stand_test_mat = standardize_test(test_mat,mean_a,std_a);
% model = svmtrain(train_label,stand_train_mat,'-t 0');
% [predicted_label, accuracy, prob_estimates] = svmpredict(test_label, stand_test_mat, model);
% res_svm_test = Evaluate(test_label,predicted_label);
% [predicted_label, accuracy, prob_estimates] = svmpredict(train_label, stand_train_mat, model);
% res_svm_train = Evaluate(train_label,predicted_label);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %projected by yang yang
% load('project_single.txt');
% samp1 = project_single(:,1:102);
% samp2 = project_single(:,103:131);
% 
% tic;
% [ new_samp1 new_samp2 ] = perm_genes( samp1,samp2,2 ); %注意实际应该是3，不是2
% toc;
% whos
% 
% tic;
% [ pvalues ] = perm_samples( new_samp1,new_samp2,5000,1); %好像1000太小，要10e9才行...内存果断爆
% toc;
% whos
% 
% tic;
% [ acc_train, acc_test ] = compute_acc( samp1, samp2, iter );
% toc;

