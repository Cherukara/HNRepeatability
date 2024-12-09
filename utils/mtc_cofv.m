function stat_cv = mtc_cofv(arr_x)
%MTC_COFV Calculate coefficient of variation of data input
%   Takes input array ARR_X, which should have N rows (where N is the number of
%   subjects) and K columns (where K is the number of repeats)
%
%   Returns a scalar STAT_CV corresponding to the within-class coefficient of
%   variation of the sample
%
%   Created by MT Cherukara, 2024-09-09

% Extract sizes
[n_subs, n_reps] = size(arr_x);

% Calculate mean across repeats and repmat it
arr_mean = repmat(mean(arr_x,'all','omitnan'),[n_subs,n_reps]);

% Calculate mean across repeats and repmat it
arr_Si = repmat(mean(arr_x,2,'omitnan'),[1,n_reps]);

% Calculate mean across subjects and repmat it
arr_Mj = repmat(mean(arr_x,1,'omitnan'),[n_subs,1]);

% Calculate sum of squares, total
arr_sst = sum((arr_x - arr_mean).^2,'all');

% Calculate sum of squares between subjects
arr_ssbs = sum((arr_mean - arr_Si).^2,'all');

% Calculate sum of squares between repeats
arr_ssbm = sum((arr_mean - arr_Mj).^2,'all');

% Calculate mean square error
stat_mse = (arr_sst - arr_ssbs - arr_ssbm)./((n_subs-1).*(n_reps-1));

% % Calculate coefficient of variation
% stat_cv = sqrt(stat_mse)./abs(arr_mean(1));

% Apply a bias correction for very small sample size
stat_cv = (1+(1/(4.*n_subs))).*sqrt(stat_mse)./abs(arr_mean(1));