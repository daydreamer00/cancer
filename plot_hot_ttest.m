function [ qvalues  ] = plot_hot_ttest(pvalues,acc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n = size(pvalues,1);
    qvalues = zeros(n,2);
    for i=1:2
%         figure(i);
        res = savePQ(pvalues(:,i),'tmp.txt',[1:n]');
%         scatter3(res(:,2),res(:,3),acc,'.');
        qvalues(:,i) = res(:,end);
    end        
    
end

