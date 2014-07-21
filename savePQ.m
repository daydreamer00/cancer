function [  ] = savePQ( pvalues,file_name,enum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [~, ~, q_values]=fdr_bh(pvalues,0.01);
%     [sorted_p sort_i] = sort(pvalues);
%     res_table = [enum sorted_p(1:n) q_values(sort_i(1:n))];
    res_table = [enum pvalues q_values];
    save(file_name,'res_table','-ascii');

end

