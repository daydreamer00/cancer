function [  ] = run_hotelling_perm_gene( k_neigh_range,proj_type_range,mean_type_range,n_perm )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    for k_neigh = k_neigh_range
        for proj_type = proj_type_range
            for mean_type = mean_type_range
                [ new_samp_pos, new_samp_neg ] = first_proj( k_neigh,proj_type,mean_type, [] );
                proj_str = get_proj_str(proj_type);
                file_name = sprintf('data/%d nn %s hotelling perm %d.txt',k_neigh,proj_str,n_perm);
                [ hot_pvalues,enum ] = hotelling_t2_test_perm_genes( new_samp_pos,new_samp_neg,n_perm);
                hot_pvalues = [enum hot_pvalues];
                save(file_name,'hot_pvalues','-ascii');
            end
        end
    end
    
    loadOriData
    
    samp1_pos = samp1(:,1:n_pos);
    samp1_neg = samp1(:,n_pos+1:end);
    samp2_pos = samp2(:,1:n_pos);% - samp1_pos;
    samp2_neg = samp2(:,n_pos+1:end);% -samp1_neg;
    
    [ hot_pvalues,enum ] = hotelling_t2_test_perm_genes( samp2_pos,samp2_neg,n_perm);
    hot_pvalues = [enum hot_pvalues];
    file_name = sprintf('data/canceronly hotelling perm %d.txt',n_perm);
    save(file_name,'hot_pvalues','-ascii');
    
    samp2_pos = samp2(:,1:n_pos) - samp1_pos;
    samp2_neg = samp2(:,n_pos+1:end) -samp1_neg;
    
    [ hot_pvalues,enum ] = hotelling_t2_test_perm_genes( samp2_pos,samp2_neg,n_perm);
    hot_pvalues = [enum hot_pvalues];
    file_name = sprintf('data/cancerminus hotelling perm %d.txt',n_perm);
    save(file_name,'hot_pvalues','-ascii');
end

