function [ new_samp_pos, new_samp_neg ,w ,thres] = k_near_proj( samp1,samp2,k_neigh,projtype,disttype,meantype,n_pos,biased,plot_index)
%K_NEAR_PROJ k neighbors project
%   k_neigh:    number of neighbors
%   projtype:   1 fisher
%               2 pSVM
%               3 libSVM
%               4 libsvm with c tuning
%   disttype:   1 euclidian distance
%   meantype:   1 mean
%               2 weighted mean
%   n_pos:      number of positive samples
%   biased:     1 yes
%               0 no
%   plot_index: genes needed to plotted
    k_neigh = k_neigh+1; %include the original one
    n1 = size(samp1,1);
    m1 = size(samp1,2);
    n2 = size(samp2,1);
    m2 = size(samp2,2);
    if(n1~=n2)
        error('not same number of genes');
    end
    
    if(m1~=m2)
        error('not same number of samples');
    end
    
    n_neg = m1 - n_pos;
    new_samp_pos = zeros(n1,n_pos);
    new_samp_neg = zeros(n1,n_neg);
    pd_squareform = zeros(n1,n1);
    w = zeros(n1,2);
    thres = zeros(n1,1);
    
    samp = [samp1 samp2];
    if(disttype == 1)
        % Euclidean dist        
        pair_dist = pdist(samp);
        pd_squareform = squareform(pair_dist);
    else
        error('invalid disttype');
    end
    
    for i = 1:n1
        [sorted,sort_i] = sort(pd_squareform(i,:));
        self_samp_pos = [samp1(i,1:n_pos); samp2(i,1:n_pos)];
        self_samp_neg = [samp1(i,n_pos+1:end); samp2(i,n_pos+1:end)];
        %Mean
        delta = (mean(pd_squareform(:))/2)^2;
        neighbor_samp_pos = zeros(n_pos,2);
        neighbor_samp_neg = zeros(n_neg,2);
        weight_sum = 0;
        for j = 1:k_neigh
            if(meantype ==1)
                weight = 1;
            elseif(meantype == 2)
                weight = exp(-pd_squareform(i,sort_i(j))^2/delta);
            else
                error('invalid meantype');
            end
            neighbor_samp_pos(:,:) = neighbor_samp_pos+weight* [samp1(sort_i(j),1:n_pos)' samp2(sort_i(j),1:n_pos)'] ;
            neighbor_samp_neg(:,:) = neighbor_samp_neg+weight* [samp1(sort_i(j),n_pos+1:end)' samp2(sort_i(j),n_pos+1:end)'];
            weight_sum = weight_sum + weight;
        end
        neighbor_samp_pos = neighbor_samp_pos/weight_sum;
        neighbor_samp_neg = neighbor_samp_neg/weight_sum;
        if(projtype == 1)
            w(i,:) = fish_proj_vec(neighbor_samp_pos,neighbor_samp_neg)';
        elseif(projtype ==2)
            w(i,:) = svm_proj(neighbor_samp_pos,neighbor_samp_neg,1)';
        elseif(projtype == 3)
            [wtmp, b] = svm_proj(neighbor_samp_pos,neighbor_samp_neg,2);
            w(i,:) = wtmp';
        elseif(projtype == 4)
            w(i,:) = svm_proj(neighbor_samp_pos,neighbor_samp_neg,3)';
        else
            error('invalid projtype');
        end
        
        if(projtype ==3 )
            thres(i) = b;
        else
            thres(i) = - w(i,:)*(mean(neighbor_samp_pos)+mean(neighbor_samp_neg))'/2;
        end
       
        if(biased)
            new_samp_pos(i,:) = w(i,:)*self_samp_pos;
            new_samp_neg(i,:) = w(i,:)*self_samp_neg;
        else
            new_samp_pos(i,:) = w(i,:)*self_samp_pos + repmat(thres(i),1,n_pos);
            new_samp_neg(i,:) = w(i,:)*self_samp_neg + repmat(thres(i),1,n_neg);
        end
        
        if (ismember(i,plot_index))
            self_samp_pos = self_samp_pos';
            self_samp_neg = self_samp_neg';
            figure(1);
%             scatter(neighbor_samp_pos(1:n_pos,1),neighbor_samp_pos(1:n_pos,2),'+','r');  
            hold on;
%             scatter(neighbor_samp_neg(1:n_neg,1),neighbor_samp_neg(1:n_neg,2),'o','r');
            scatter(self_samp_pos(:,1),self_samp_pos(:,2),'+','r');
            scatter(self_samp_neg(:,1),self_samp_neg(:,2),'o','k','fill');
            min_x = min([self_samp_pos(1:n_pos,1); self_samp_neg(1:n_neg,1)]);
            max_x = max([self_samp_pos(1:n_pos,1); self_samp_neg(1:n_neg,1)]);
            centerx = (mean(neighbor_samp_pos)+mean(neighbor_samp_neg))'/2;
%             plotbound(centerx,w(i,:)');
            plotbound2(thres(i),w(i,:)');
            hold off;
            
            figure(2);
            plotHist2(new_samp_pos(i,:),new_samp_neg(i,:));
            waitforbuttonpress;
        end
        
    end
           
end

function [] = plotbound(x,w)
    wb = [w(2,1); -w(1,1)];
    x2 = x+wb;
    x1 = x-wb;
    line([x1(1,1) x2(1,1)],[x1(2,1),x2(2,1)]);
end

function [] = plotbound2(t,w)
    xrang = get(gca,'xlim');
    yrang = get(gca,'ylim');
    wb = [w(2,1); -w(1,1)];
    m = w(1,1)/(-w(2,1));
    b = t/(-w(2,1));
    x = xrang(1):0.00001:xrang(2);
    y = m*x+b;
    index = x>=xrang(1) & x<xrang(2) & y < yrang(2) & y >=xrang(1);
    plot(x(index),y(index));
%     plot(x,y);
end

