function [w,model] = svmtrain_prj(labels,samp,svmoption)
    bestcv = 0;
    c_range = 1:10;
    g_range = 1;
    [n m] = size(samp);
    %             error_his = zeros(length(c_range),length(g_range));
    for log2c = c_range,
        for log2g = g_range
            cmd = [svmoption ' -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            rand_i = randperm(n);
            rand_samp = samp(rand_i,:);
            n_fold = 5;
            iter = 100;
            cv = 0;
            
            for j = 1:iter
                cv_i = crossvalind('Kfold',labels,n_fold);
                
                for i = 1:n_fold
                    test_i = (cv_i == i); train_i = ~test_i;
                    model = svmtrain(labels(train_i), samp(train_i,:), cmd);
                    [predicted_label, accuracy, prob_estimates] = svmpredict(labels(test_i), samp(test_i,:), model,'-q');
                    res_svm = Evaluate(labels(test_i),predicted_label);
                    mean_acc = mean(res_svm(2:3));
                    cv = cv+ mean_acc;
                end
                
            end
            
            cv = cv/n_fold/iter;
            if (cv > bestcv),
                bestcv = cv; bestlog2c = log2c; bestlog2g = log2g;
                best_cmd = [svmoption ' -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            end
            
            i = log2c-c_range(1)+1;
            j = log2g-g_range(1)+1;
            %                     error_his(i,j) = cv;
        end
    end
    svmoption = best_cmd;
    model = svmtrain(labels(train_i), samp(train_i,:), svmoption);
    w = model.SVs' * model.sv_coef;
    fprintf('(best log2c=%g, log2g=%g, best cv=%g)\n', bestlog2c, bestlog2g, bestcv);
end