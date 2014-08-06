function [scaled] = scale_data(data)
    maxv = max(data);
    minv = min(data);
    scaled = (data - minv)/(maxv-minv);
end