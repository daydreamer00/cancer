function [ pvalues ] = getPFromT( t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n = size(t,1);
    pvalues = zeros(n,1);
    for i=1:n
        t_statistic = t(i,1);
        if t_statistic~=0
            df = (V1+V2)^2/(V1^2/(n1_rev-1)+V2^2/(n2_rev-1)+eps);
        else
            df = 1;
        end
        p_value(i,1) = 1-tcdf(t_statistic,df);
        p_value(i,1) = p_value(i,1)*2;   % Bilateral distribution
    end

end

