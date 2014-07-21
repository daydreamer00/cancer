function [ error_his ] = gridsearch()
%GRIDSEARCH 
%  Usage: error_his = gridsearch()

    ori_data = xlsread('miRNA_with_metastasis.xlsx');
    samp1 = ori_data(5:end,ori_data(1,:)==1);
    samp2 = ori_data(5:end,ori_data(1,:)==2);
    n_pos = sum(ori_data(4,:)==0)/2;
    n_neg = size(samp1,2) - n_pos;
    train_r = 0.7;
    n_train_pos = fix(n_pos*train_r);
    n_train_neg = fix(n_neg*train_r);
    train_label = [ ones(n_train_pos,1); zeros(n_train_neg,1)];
    n_test_pos = n_pos - n_train_pos;
    n_test_neg = n_neg - n_train_neg;
    test_label = [ ones(n_test_pos,1); zeros(n_test_neg,1)];
    samp = samp1';
    train_mat = samp([1:n_train_pos n_pos+(1:n_train_neg)],:);
    test_mat = samp([n_train_pos+(1:n_test_pos) n_pos+n_train_neg+(1:n_test_neg)],:);
    [stand_train_mat mean_a std_a ]  = standardize(train_mat);
    stand_test_mat = standardize_test(test_mat,mean_a,std_a);
    
    bestcv = 0;
    c_range = 0:10;
    g_range = -15:-4;
    error_his = zeros(length(c_range),length(g_range));
    for log2c = c_range,
        for log2g = g_range
            cmd = ['-q -t 2 -v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            cv = svmtrain(train_label, stand_train_mat, cmd);
            if (cv > bestcv),
                bestcv = cv; bestlog2c = log2c; bestlog2g = log2g;
                best_cmd = ['-q -t 2 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            end
            i = log2c-c_range(1)+1;
            j = log2g-g_range(1)+1;
            error_his(i,j) = cv;
        end
    end
    fprintf('(best log2c=%g, log2g=%g, err_rate=%g)\n', bestlog2c, bestlog2g, bestcv);
    surf(g_range,c_range,error_his);
    
    model = svmtrain(train_label,stand_train_mat,best_cmd);
    [predicted_label, accuracy, prob_estimates] = svmpredict(test_label, stand_test_mat, model);
    res_svm_test = Evaluate(test_label,predicted_label);
end

