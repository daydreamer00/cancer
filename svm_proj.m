function [ w , b] = svm_proj( samp1,samp2,svm_type )
%svm_proj Summary of this function goes here
%   svm_type:   1 psvm
%               2 libsvm
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(m1~=m2)
        error('not same number of dim');
    end
    
    labels = ones(n1+n2,1);
    labels(n1+1:end,1) = -labels(n1+1:end,1);
    samp = [samp1;samp2];
    mean_a = [0 0]; std_a = [1 1];
%     [ samp ,mean_a, std_a] = standardize( samp );
    switch svm_type
        case 1
            w = psvm(samp,labels,0,1);
            w = w./std_a';
        case 2
            [ samp ,mean_a, std_a] = standardize( samp );
%             bestcv = 0;
%             c_range = 10:20;
%             g_range = 1;
% %             error_his = zeros(length(c_range),length(g_range));
%             for log2c = c_range,
%                 for log2g = g_range
%                     cmd = ['-s 0 -q -t 0 -v 10 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
%                     cv = svmtrain(labels, samp, cmd);
%                     if (cv > bestcv),
%                         bestcv = cv; bestlog2c = log2c; bestlog2g = log2g;
%                         best_cmd = ['-s 0 -q -t 2 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
%                     end
%                     i = log2c-c_range(1)+1;
%                     j = log2g-g_range(1)+1;
% %                     error_his(i,j) = cv;
%                 end
%             end
%             svmoption = best_cmd;
%             fprintf('(best log2c=%g, log2g=%g, err_rate=%g)\n', bestlog2c, bestlog2g, bestcv);
            svmoption = sprintf('-s 0 -t 0 -w1 %f -q',n2/n1);
%             svmoption = sprintf('-s 0 -t 0 -w1 %f -q',1);
            model = svmtrain(labels,samp,svmoption);
            
            w = model.SVs' * model.sv_coef;
            
            b = (-model.rho-w'*(mean_a./std_a)');
            w = w./std_a';
            
%             [predicted_label, accuracy, prob_estimates] = svmpredict(labels, samp, model);
%             res_svm = Evaluate(labels,predicted_label);
%             dec = samp*w+b;
            b = b/sqrt(w'*w);
        case 3
            [ samp ,mean_a, std_a] = standardize( samp );
            svmoption = sprintf('-s 0 -t 0 -w1 %f -q',n2/n1);
%             svmoption = sprintf('-s 0 -t 0 -w1 %f -q',1);
            w = svmtrain_prj(labels,samp,svmoption);
%             b = (-model.rho-w'*(mean_a./std_a)');
            w = w./std_a';
        otherwise
            error('wrong svm type');
    end
    
    w = w / sqrt(w'*w);    

end

