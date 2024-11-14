function t_anova = mtc_anova_subs(arr_x)
%MTC_ANOVA_SUBS Perform repeated-measures ANOVA to test for subject-level
%differences
%
%   Takes input array ARR_X, which should have N rows (where N is the number of
%   subjects) and K columns (where K is the number of repeats)
%
%   Returns a table T_ANOVA which is the ANOVA results
%
%   Created by MT Cherukara, 2024-09-09

% Extract sizes
n_reps = size(arr_x,2);

% Create an array of subject IDs (this will be a dummy variable in the ANOVA)
arr_sub = repelem([{'a'},{'b'},{'c'},{'d'},{'e'},{'f'},{'g'},{'h'},{'i'},{'j'}],n_reps)';

% Create a table
t_x = table(arr_sub(:),arr_x(:),'VariableNames',{'subs','measure'});
t_mes = table(1,'VariableNames',{'Measurements'});

% Create the repeated measures model
m_rm = fitrm(t_x,'measure~subs','WithinDesign',t_mes);

% Do ANOVA
t_anova = anova(m_rm);