function [  ] = get_smallest_pq(k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dir_name = 'data/pq/';
    dir_files = dir(dir_name);
    n_files = size(dir_files,1);
    for i = 1:n_files
        i
        file_name = dir_files(i).name;
        if(regexp(file_name,'\.txt\s*$'))
            p_data = load([dir_name file_name]);
            [m, n] = size(p_data);
            sorted = sortrows(p_data,[n n-1]);
            sorted = sorted(1:k,:);
            save([dir_name 'smallest/' file_name],'sorted','-ascii');
        end
    end

end

