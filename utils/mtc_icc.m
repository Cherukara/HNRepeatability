function stat_icc = mtc_icc(arr_x)
%MTC_ICC Calculate intra-class correlation coefficient of data
%   Takes input array ARR_X, which should have N rows (where N is the number of
%   subjects) and K columns (where K is the number of repeats)
%
%   Returns a scalar STAT_ICC corresponding to the intra-class correlation
%   coefficient (one-way model)
%
%   Created by MT Cherukara, 2024-09-09

% Extract sizes
[n_subs, n_reps] = size(arr_x);

% Calculate mean across repeats and repmat it
arr_mean = repmat(mean(arr_x,'all','omitnan'),[n_subs,n_reps]);

% Calculate mean across repeats and repmat it
arr_Si = repmat(mean(arr_x,2,'omitnan'),[1,n_reps]);

% Calculate mean squares within subjects
arr_msws = sum((arr_Si - arr_x).^2,'all')./(n_subs.*(n_reps-1));

% Calculate mean squares between subjects
arr_msbs = sum((arr_Si - arr_mean).^2,'all')./(n_subs-1);

% Calculate ICC
stat_icc = (arr_msbs - arr_msws) ./ (arr_msbs + (n_reps-1).*arr_msws);

