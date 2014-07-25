function [ new_samp_pos, new_samp_neg ] = first_proj( k_neigh_range,proj_type_range,mean_type_range, plot_indexs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   proj_type:  1:fisher 2:pSVM 3:libSVM 4: tuned libSVM
%   mean_type:  1:normal 2:weighted
    loadOriData   
    
    for k_neigh = k_neigh_range
        for proj_type = proj_type_range
            for mean_type = mean_type_range
                [ new_samp_pos, new_samp_neg ] = k_near_proj( samp1,samp2,k_neigh,proj_type,1,mean_type,n_pos,1,plot_indexs);
                new_samp = [new_samp_pos new_samp_neg];
                if(proj_type ==1 )
                    save(sprintf('%d nn fisher.txt',k_neigh),'new_samp','-ascii');
                elseif(proj_type==2)
                    save(sprintf('%d nn psvm.txt',k_neigh),'new_samp','-ascii');
                elseif(proj_type==3)
                    save(sprintf('%d nn libsvm.txt',k_neigh),'new_samp','-ascii');
                end
            end
        end
    end  

end

