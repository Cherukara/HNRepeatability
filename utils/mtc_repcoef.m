function stat_rc = mtc_repcoef(arr_x)
%MTC_REPCOEF Calculate the repeatability coefficient of data input
%   Takes input array ARR_X, which should have N rows (where N is the number of
%   subjects) and K columns (where K is the number of repeats)
%
%   Returns a scalar STAT_RC corresponding to the sample repeatability
%   coefficient, based on a 95% confidence interval
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
arr_ssbs = sum((arr_Si - arr_mean).^2,'all');

% Calculate sum of squares between repeats
arr_ssbm = sum((arr_Mj - arr_mean).^2,'all');

% Calculate mean square error
stat_mse = (arr_sst - arr_ssbs - arr_ssbm)./((n_subs-1).*(n_reps-1));

% Calculate repeatability coefficient
stat_rc = 1.96.*sqrt(2.*stat_mse);

end