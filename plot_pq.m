function [ output_args ] = plot_pq( k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    dir_name = 'data/pq/smallest/';
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
            scatter(pq_data(1:k,n-1),pq_data(1:k,n),'.');
            title(file_name)
            i_fig = i_fig+1;
        end
    end

end

