function [ output_args ] = rank_first_prj(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    sheet_count = 0;
    
    for proj_type = [1 3]
        for k_neigh = 0:3
            for mean_type = 1
                prj_str = get_proj_str(proj_type);
                sheet_name = sprintf('%dnn %s',k_neigh,prj_str);
                table = xlsread('first_prj.xls',sheet_name);
                i_t_p = 3;
                i_t_q = 5;
                i_acc = 8;
                sorted = sortrows(table,i_t_p);
                
                n = size(table,1);
                ratio_top = 0.1;
                n_top = fix(n*ratio_top);
                                
                err_rate = 1- sorted(1:n_top,i_acc);
                err_rate = scale_data(err_rate);
                pvalue = scale_data(sorted(1:n_top,i_t_p));
                qvalue = scale_data(sorted(1:n_top,i_t_q));
                d_sq = err_rate.*err_rate+pvalue.*pvalue+qvalue.*qvalue;
                id_arr = sorted(1:n_top,1);
                rank_table = [id_arr d_sq];
                rank_table = sortrows(rank_table,2);
                
%                 scatter3(pvalue,qvalue,err_rate);
%                 waitforbuttonpress;
                
                sheet_count = sheet_count+1;
                
                file_name = 'first_prj_rank.xls';
                xlswrite(file_name,rank_table,sheet_name);
                start_col = char((sheet_count-1)*2+'A');
                xlswrite(file_name,rank_table(1:10,:),'top10',[start_col '2']);
            end

        end
    end
end

function [scaled] = scale_data(data)
    maxv = max(data);
    minv = min(data);
    scaled = (data - minv)/(maxv-minv);
end
