
function [h crit_p adj_p]=fdr_bh_no_min(pvals,q,m_all)

% if nargin<1,
%     error('You need to provide a vector or matrix of p-values.');
% else
%     if ~isempty(find(pvals<0,1)),
%         error('Some p-values are less than 0.');
%     elseif ~isempty(find(pvals>1,1)),
%         error('Some p-values are greater than 1.');
%     end
% end

if nargin<2,
    q=.05;
end

s=size(pvals);
if (length(s)>2) || s(1)>1,
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[~, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests
if nargin == 2
    m_all = m;
end

%BH procedure for independence or positive dependence
thresh=(1:m)*q/m_all;
wtd_p=m_all*p_sorted./(1:m);

if nargout>2,
    %compute adjusted p-values
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    index = find(wtd_p_sorted>1);
    if size(index,2) == size(wtd_p_sorted,2)
        wtd_p_sorted = ones(size(wtd_p_sorted));
    elseif 1~=isempty(index)
        wtd_p_sorted(index) = wtd_p_sorted(index(1)-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    adj_p = wtd_p;
    adj_p=reshape(adj_p(unsort_ids),s);
end

rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
end