function [ string ] = get_proj_str( proj_type )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    switch proj_type
        case 1
            string = 'fisher';
        case 2
            string = 'psvm';
        case 3
            string = 'libsvm';
        otherwise
            error('wrong proj_type');
    end

end

