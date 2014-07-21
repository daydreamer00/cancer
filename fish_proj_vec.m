function [w] = fish_proj_vec(mat1,mat2)
    n1 = size(mat1,1);
    n2 = size(mat2,1);
    mu1 = mean(mat1);
    mu2 = mean(mat2);
    B_mat = (mu1 - mu2)'*(mu1-mu2);
%     cov1 = (mat1-repmat(mu1,n1,1))'*(mat1-repmat(mu1,n1,1));
%     cov2 = (mat2-repmat(mu2,n2,1))'*(mat2-repmat(mu2,n2,1));
    cov1 = cov(mat1);
    cov2 = cov(mat2);
    W_mat = (cov1*(n1-1)+cov2*(n2-1));
    w = (mu1-mu2)/W_mat;
    w = w' / sqrt(w*w');
end