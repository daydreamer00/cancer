function [ stand_mat mean_a std_a] = standardize( data_mat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    m = size(data_mat,1);
    mean_a= mean(data_mat,1);
    mean_mat = ones(m,1)*mean_a;
    std_a = std(data_mat,0,1);
    std_mat = ones(m,1)*std_a;
    stand_mat = (data_mat - mean_mat)./std_mat;
end