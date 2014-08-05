function [ output_args ] = run_perm_genes_compute_acc( k_neigh_range,proj_type_range,mean_type_range,n_perm,iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    for k_neigh = k_neigh_range
        for proj_type = proj_type_range
            for mean_type = mean_type_range
                [ new_samp_pos, new_samp_neg ] = first_proj( k_neigh,proj_type,mean_type, [] );
                proj_str = get_proj_str(proj_type);                
                [ res_table ] = compute_acc_perm_genes( new_samp_pos,new_samp_neg,n_perm,iter);
                
                file_name = '2comb_acc.xls';
                prj_str = get_proj_str(proj_type);
                header = [repmat({'gene id'},[n_perm 1]),{'acc train','acc test','avg acc'}];
                xlsdata = [header;num2cell(res_table)];                
                sheet_name = sprintf('%dnn %s',k_neigh,prj_str);
                
                xlswrite(file_name,xlsdata,sheet_name);
            end
        end
    end
    
    loadOriData
    
    samp1_pos = samp1(:,1:n_pos);
    samp1_neg = samp1(:,n_pos+1:end);
    samp2_pos = samp2(:,1:n_pos);% - samp1_pos;
    samp2_neg = samp2(:,n_pos+1:end);% -samp1_neg;
    
     [ res_table ] = compute_acc_perm_genes( samp2_pos,samp2_neg,n_perm,iter);
    
    prj_str = 'cancer_only';
     header = [repmat({'gene id'},[n_perm 1]),{'acc train','acc test','avg acc'}];
     xlsdata = [header;num2cell(res_table)];
     sheet_name = sprintf('%s',prj_str);
     xlswrite(file_name,xlsdata,sheet_name);
    
    
    samp2_pos = samp2(:,1:n_pos) - samp1_pos;
    samp2_neg = samp2(:,n_pos+1:end) -samp1_neg;
    
     [ res_table ] = compute_acc_perm_genes( samp2_pos,samp2_neg,n_perm,iter);
    
    prj_str = 'cancer_minus';
    header = [repmat({'gene id'},[n_perm 1]),{'acc train','acc test','avg acc'}];
     xlsdata = [header;num2cell(res_table)];
     sheet_name = sprintf('%s',prj_str);
     xlswrite(file_name,xlsdata,sheet_name);

end

