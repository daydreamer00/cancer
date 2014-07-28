function [ output_args ] = run_compute_first_proj(k_neigh_range,proj_type_range,mean_type_range  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    loadOriData
    for k_neigh = k_neigh_range
        for proj_type = proj_type_range
            for mean_type = mean_type_range
                [ new_samp_pos, new_samp_neg ] = k_near_proj( samp1,samp2,k_neigh,proj_type,1,mean_type,n_pos,1,[]);
                new_samp = [new_samp_pos new_samp_neg];
                
                [ ttest_pvalues ] = welch_t_test( new_samp_pos,new_samp_neg );
                
                [ new_samp1, new_samp2, enum] = perm_genes(  new_samp_pos,new_samp_neg,1 );
                
                [ acc_train, acc_test ,train_eval, test_eval] = compute_acc_sep( samp1,samp2,k_neigh,proj_type,1,mean_type,n_pos,1000 );
                avg_acc = (acc_train+acc_test)/2;
                acc_sep = [acc_train acc_test avg_acc];
                [ acc_train, acc_test ] = compute_acc( new_samp1,new_samp2,500 );
                avg_acc = (acc_train+acc_test)/2;
                acc = [acc_train acc_test avg_acc];
                
                %hotelling
                [ hot_pvalues ] = hotelling_t2_test_batch( samp1,samp2,n_pos );
                
                qvalues = plot_hot_ttest([hot_pvalues ttest_pvalues],ones(n1,1)-acc_sep(:,2));
                
                res_table_perm = [[1:n1]' hot_pvalues ttest_pvalues qvalues acc_sep train_eval test_eval acc];
                
                file_name = 'first_prj.xls';
                prj_str = get_proj_str(proj_type);
                header = {'gene id','hot pvalue','ttest pvalue','hot qvalue','ttest qvalue','acc train','acc test','avg acc','accuracy','sensitivity','specificity','precision','recall','f_measure','gmean','accuracy','sensitivity','specificity','precision','recall','f_measure','gmean','','',''};
                xlsdata = [header;num2cell(res_table_perm)];                
                sheet_name = sprintf('%dnn %s',k_neigh,prj_str);
                
                xlswrite(file_name,xlsdata,sheet_name);
                
            end
        end
    end  
    

end

