function [ output_args ] = run_compute_q(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    dir_name = 'data/';
    dir_files = dir(dir_name);
    n_files = size(dir_files,1);
    for i = 1:n_files
        i
        file_name = dir_files(i).name;
        if(regexp(file_name,'\.txt\s*$'))
            p_data = load([dir_name file_name]);
            [m, n] = size(p_data);
            enum = p_data(:,1:n-1);
            pvalues = p_data(:,n);
            [~, ~, qvalues]=fdr_bh(pvalues,0.01);
            pqdata = [p_data qvalues];
            save([dir_name 'pq/' file_name],'pqdata','-ascii');
        end
    end
end

