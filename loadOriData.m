%load data
ori_data = xlsread('miRNA_with_metastasis.xlsx');
save('ori_data.mat','ori_data');
load('ori_data.mat');
samp1 = ori_data(5:end,ori_data(1,:)==1);
samp2 = ori_data(5:end,ori_data(1,:)==2);

n1 = size(samp1,1);

n_pos = sum(ori_data(4,:)==0)/2;