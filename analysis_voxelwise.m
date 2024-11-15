% ANALYSIS_VOXELWISE.m calculates voxel-wise repeatability metrics for QSM
%
% MT Cherukara, October 2024

clearvars;
close all;

% Data directory
dir_data = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/HNRepeatability_Data/';

% Methods we are comparing
meth_names = {'LBV','PDF','VSHARP'};
% meth_names = {'LPU','SEGUE'};
n_meth = length(meth_names);

% Long ROI names
names_roi = {'Thalamus','Caudate Nucleus','Putamen','Globus Pallidus'};
n_rois = length(names_roi);

% Lengths
n_subs = 10;
n_reps = 6;

% We might want to rename the methods
names_methnice = meth_names;

% Pre-allocate storage arrays, averaging over subjects and sessions, storing
% mean and standard deviation
stat_rmse_roi = zeros(n_meth,n_rois,2);
stat_xsim_roi = zeros(n_meth,n_rois,2);

stat_rmse_all = zeros(n_meth,n_subs,n_reps-1,n_rois);
stat_xsim_all = zeros(n_meth,n_subs,n_reps-1,n_rois);


%% Loop through methods, load and fill in the table
for mm = 1:n_meth

    % Load the Data
    struct_m = load(strcat(dir_data,'Voxel_Repeatability_',meth_names{mm},'_data.mat'));

    % Fill in the mean for each one
    stat_rmse_roi(mm,:,1) = mean(struct_m.res_rmse(:,2:n_reps,1:n_rois),[1,2]);
    stat_xsim_roi(mm,:,1) = mean(struct_m.res_xsim(:,2:n_reps,1:n_rois),[1,2]);

    % Fill in the standard deviation for each one
    stat_rmse_roi(mm,:,2) = std(struct_m.res_rmse(:,2:n_reps,1:n_rois),[],[1,2]);
    stat_xsim_roi(mm,:,2) = std(struct_m.res_xsim(:,2:n_reps,1:n_rois),[],[1,2]);

    % Put all the data in a big array so that we can do Repeated-Measures ANOVA
    %       Dimensions: METHODS * SUBJECTS * REPEATS * ROIs
    stat_rmse_all(mm,:,:,:) = struct_m.res_rmse(:,2:n_reps,1:n_rois);
    stat_xsim_all(mm,:,:,:) = struct_m.res_xsim(:,2:n_reps,1:n_rois);

end % for mm = 1:n_meth


%% One-way ANOVA

% Do one-way ANOVA
res_anova = anova(stat_xsim_roi(:,:,1)');

% Pull p-value
res_pval = res_anova.stats.pValue(1);

% If there is a significance, do multiple comparisons
if res_pval < 0.05

    t_mult_r = multcompare(res_anova);

    % Pull out only the statistically significant pairs
    t_mult_r = t_mult_r(t_mult_r.pValue < 0.05, :);

    % Loop over remaining pairs
    for pp = 1:size(t_mult_r,1)
        disp([names_methnice{t_mult_r.Group1(pp)},' - ',names_methnice{t_mult_r.Group2(pp)},...
             ' (ROI) are significantly different, at p = ',num2str(t_mult_r.pValue(pp))]);
    end


else % if res_pval < 0.05
    disp('No significant differences (RMSE)');

end % if res_pval < 0.05 // else


%% Repeated Measures ANOVA

% We need a vector of "subject names" from A (1) to J (10) that repeats 5 times
arr_sub = repelem([{'a'},{'b'},{'c'},{'d'},{'e'},{'f'},{'g'},{'h'},{'i'},{'j'}],n_reps-1)';

% Pre-allocate an array for storing ROI p-Values
res_ranova = zeros(n_rois,1);

% Numbers for pairwise comparisons
pairs_mult = nchoosek(1:n_meth,2);
n_pairs = size(pairs_mult,1);

% Pre-allocate a table for pairwise comparisons
t_pairs = table('Size',[n_pairs,2+n_rois],...
                'VariableTypes',[{'string','string'},repelem({'logical'},n_rois)],...
                'VariableNames',[{'Method1','Method2'},names_roi]);

t_pairs.Method1 = [names_methnice(pairs_mult(:,1))]';
t_pairs.Method2 = [names_methnice(pairs_mult(:,2))]';

% Loop over ROIs and do the ANOVA in each one
for rr = 1:n_rois

    % Pull out the current ROI
    arr_stat_roi = stat_rmse_all(:,:,:,rr);

    % Create table for this ROI
    t_roi = cell2table(arr_sub,'VariableNames',{'sub'});

    % Loop over methods and add columns
    for mm = 1:n_meth

        % Pull out the current method
        arr_stat_meth = squeeze(arr_stat_roi(mm,:,:))';

        % Add to the table
        t_roi = addvars(t_roi,arr_stat_meth(:),'NewVariableNames',{sprintf('meas%d',mm)});

    end % for mm = 1:n_meth

    % Create measurements table
    t_meas = table((1:n_meth)','VariableNames',{'Measurements'});

    % Create Repeated Measures Model
    rm = fitrm(t_roi,'meas1-meas3~sub','WithinDesign',t_meas);

    % Now, do the repeated measures ANOVA
    t_ranova = ranova(rm);

    % We want the lower-bound p-Value from the Intercept:Measurements term,
    % because this tells us whether the "Measurements" (methods) are
    % significantly different
    res_ranova(rr) = t_ranova.pValueLB(1);

    % Pairwise comparison
    t_mult = multcompare(rm,'Measurements');

    % Loop through pairs and fill in the pairwise comparison table
    for pp = 1:n_pairs

        % Find the correct row
        i_row = find(t_mult.Measurements_1 == pairs_mult(pp,1) & ...
                     t_mult.Measurements_2 == pairs_mult(pp,2) );

        % Insert p-value
        t_pairs.(names_roi{rr})(pp) = t_mult.pValue(i_row) < 0.05;

    end % for pp = 1:n_pairs

end % for rr = 1:n_rois



%% Bar Chart of RMSE

% Figure
f1 = figure(1); clf;
set(f1,'Position',[150,300,(200 + 100.*n_rois),550]);

% Bar Chart
br1 = bar(stat_rmse_roi(:,:,1)',1,'FaceColor','Flat');
box on; hold on;

% Pre-allocate array of x positions for the ends of the boxes
xpos = zeros(n_rois,n_meth);

for ee = 1:n_meth
    xpos(:,ee) = br1(ee).XEndPoints;
end

% Pull out centre points and ends of the error bars
br_c = stat_rmse_roi(:,:,1)';
br_e = stat_rmse_roi(:,:,2)';

% Plot the error bars
er1 = errorbar(xpos(:),br_c(:),br_e(:),'k','LineStyle','none','LineWidth',1);

% Labels
ylim([0,1]);
ylabel('NRMSE');
xticks(1:n_rois);
xticklabels(names_roi);
% 
% % Significance bars - whole brain
% for pp = 1:size(t_mult_w,1)
%     plot([xpos(t_mult_w.Group1(pp),1),xpos(t_mult_w.Group2(pp),1)],...
%          [110+(4*pp),110+(4*pp)],...
%          'k-','LineWidth',1.25)
% end
% 
% % Significance bars - rois
% for pp = 1:size(t_mult_r,1)
%     plot([xpos(t_mult_r.Group1(pp),2),xpos(t_mult_r.Group2(pp),2)],...
%          [122+(4*pp),122+(4*pp)],...
%          'k-','LineWidth',1.25)
% end

% Legend
legend(names_methnice,'Location','NorthEast')
legend('boxoff');


%% Bar Chart of XSIM

% Figure
f2 = figure(2); clf;
set(f2,'Position',[250,300,(200 + 100.*n_rois),550]);

% Bar Chart
br2 = bar(stat_xsim_roi(:,:,1)',1,'FaceColor','Flat');
box on; hold on;

% Pre-allocate array of x positions for the ends of the boxes
xpos = zeros(n_rois,n_meth);

for ee = 1:n_meth
    xpos(:,ee) = br2(ee).XEndPoints;
end

% Pull out centre points and ends of the error bars
br_c = stat_xsim_roi(:,:,1)';
br_e = stat_xsim_roi(:,:,2)';

% Plot the error bars
er2 = errorbar(xpos(:),br_c(:),br_e(:),'k','LineStyle','none','LineWidth',1);

% Labels
ylim([0,1]);
ylabel('XSIM');
xticks(1:n_rois);
xticklabels(names_roi);

% % Significance bars - whole brain
% for pp = 1:size(t_mult_w,1)
%     plot([xpos(t_mult_w.Group1(pp),1),xpos(t_mult_w.Group2(pp),1)],...
%          [110+(4*pp),110+(4*pp)],...
%          'k-','LineWidth',1.25)
% end
% 
% % Significance bars - rois
% for pp = 1:size(t_mult_r,1)
%     plot([xpos(t_mult_r.Group1(pp),2),xpos(t_mult_r.Group2(pp),2)],...
%          [122+(4*pp),122+(4*pp)],...
%          'k-','LineWidth',1.25)
% end

% Legend
legend(names_methnice,'Location','NorthWest')
legend('boxoff');



