function [  ] = write_to_file( mat, filename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    save(filename,'mat','-ascii');
end

