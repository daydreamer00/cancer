function [  ] = plotHist2( x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    min_v = min([x y]);
    max_v = max([x y]);
    nbin = 15;
    d = (max_v - min_v)/nbin;
    bin_x = min_v:d:max_v;
    [fx] = histc(x,bin_x);
    [fy] = histc(y,bin_x);
    Y = [fx/sum(fx); fy/sum(fy)]';
    bar(bin_x',Y,1);
%     hold on;
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','r','EdgeColor','w');
%     [f] = histc(y,bin_x);
%     bar(bin_x,f/sum(f));
%     hold off;
end

