function [ pvalues ] = perm_samples( samp1,samp2,iter_n,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    
    samp1_mean = mean(samp1(:,:),2);
    samp2_mean = mean(samp2(:,:),2);
    samp1_var = var(samp1(:,:),0,2);
    samp2_var = var(samp2(:,:),0,2);
    t = abs((samp1_mean - samp2_mean)./sqrt(samp1_var/m1+samp2_var/m2));
      
    samp = [samp1 samp2];
    
    %t_perm = single(zeros(n1,iter_n));
    n_larger_than_t = zeros(1,n1);
    
    batch_size = 10;
    n_batch = single(zeros(batch_size,n1));
      
%     for k = 1:ceil(iter_n/batch_size)
%         n_batch = single(zeros(batch_size,n1));
%         disp(sprintf('batch %d start\n',i));
%         for i = 1:batch_size
%             disp(sprintf('in batch %d start\n',i));
%             index = randperm(m1+m2);
%             samp1_mean_perm = mean(samp(:,index(1:m1)),2);
%             samp2_mean_perm = mean(samp(:,index(m1+1:m1+m2)),2);
%             samp1_var_perm = var(samp(:,index(1:m1)),0,2);
%             samp2_var_perm = var(samp(:,index(1:m1)),0,2);
%             %t_perm(:,i) = abs((samp1_mean - samp2_mean)./sqrt(samp1_var/m1+samp2_var/m2));
%             t_perm = abs((samp1_mean_perm - samp2_mean_perm)./sqrt(samp1_var_perm/m1+samp2_var_perm/m2));
%             t_perm = sort(t_perm(:));
%             n_batch(i,:) = (n1 - bsearch(t_perm,t)');
%             disp(sprintf('in batch %d end\n',i));
%         end
%         n_larger_than_t = n_larger_than_t+sum(n_batch);
%     end
    
    for i = 1:iter_n
%         fprintf('iter %d start\n',i);
        index = randperm(m1+m2);
        
        new_samp1 = samp(:,index(1:m1));
        new_samp2 = samp(:,index(m1+1:m1+m2));
        
        samp1_mean_perm = zeros(1,n1);
        samp2_mean_perm = zeros(1,n1);
        samp1_var_perm = zeros(1,n1);
        samp2_var_perm = zeros(1,n1);
%         tic;
        samp1_mean_perm = mean(samp(:,index(1:m1)),2);
        samp2_mean_perm = mean(samp(:,index(m1+1:m1+m2)),2);
        samp1_var_perm = var(samp(:,index(1:m1)),0,2);
        samp2_var_perm = var(samp(:,index(m1+1:m1+m2)),0,2);
%         toc;
%         tic;
%         parfor j = 1:n1
%             samp1_mean_perm(j) = mean(new_samp1(j,:));
%             samp2_mean_perm(j) = mean(new_samp2(j,:));
%             samp1_var_perm(j) = var(new_samp1(j,:));
%             samp2_var_perm(j) = var(new_samp2(j,:));
%         end
%         toc;
%         tic;
        %t_perm(:,i) = abs((samp1_mean - samp2_mean)./sqrt(samp1_var/m1+samp2_var/m2));
        t_perm = abs((samp1_mean_perm - samp2_mean_perm)./sqrt(samp1_var_perm/m1+samp2_var_perm/m2));
        
        if(type == 1)
            t_perm = sort(t_perm(:));
%             fprintf('iter %d start search\n',i);
            n_larger_than_t = n_larger_than_t+(n1 - bsearch(t_perm,t)');
        elseif(type == 2)
            n_larger_than_t = n_larger_than_t+(t_perm>=t)';
        else
            error('invalid type');
        end
%         n_larger_than_t = n_larger_than_t+arrayfun(@(x) sum(t_perm>x),t);
%         fprintf('iter %d end\n',i);
%         toc;
    end 
    if(type==1)
        N = iter_n*n1;
    elseif(type ==2)
        N = iter_n;
    end
%     t_perm = sort(t_perm(:));
%     pvalues = (N - bsearch(t_perm,t))/N;
    n_larger_than_t(n_larger_than_t==0) = 1;
    pvalues = n_larger_than_t'/N;
    
end

