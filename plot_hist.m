function [ output_args ] = plot_hist( k,nbin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    nbin = 400;
    dir_name = 'data3/pq/smallest/';
    dir_files = dir(dir_name);
    n_files = size(dir_files,1);
    i_fig = 1;
    for i = 1:n_files
        i
        file_name = dir_files(i).name;
        if(regexp(file_name,'\.txt\s*$'))
            display(file_name);
            pq_data = load([dir_name file_name]);
            [m, n] = size(pq_data);
            k = min([m k]);
            figure(i_fig);
            i_data = [];
            for j = 2:n-2
                i_data = [i_data; pq_data(:,j)];              
            end
            hist(i_data,nbin);
            title(file_name)
            i_fig = i_fig+1;
        end
    end

end
