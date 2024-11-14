% ANALYSIS_VOXELWISE.m calculates voxel-wise repeatability metrics for QSM
%
% MT Cherukara, October 2024

clearvars;
close all;

% Data directory
dir_data = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/HNRepeatability_Data/';

% Methods we are comparing
meth_names = {'LPU','SEGUE'};
n_meth = length(meth_names);

% Long ROI names
names_roi = {'Thalamus','Caudate Nucleus','Putamen','Globus Pallidus','Whole Brain'};
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
% gr_rmse_roi = zeros(n_meth,2);
% gr_xsim_roi = zeros(n_meth,2);


%% Loop through methods, load and fill in the table
for mm = 1:n_meth

    % Load the Data
    struct_m = load(strcat(dir_data,'Voxel_Repeatability_',meth_names{mm},'_data.mat'));

    % Fill in the mean for each one
    stat_rmse_roi(mm,:,1) = mean(struct_m.res_rmse(:,2:n_reps,:),[1,2]);
    stat_xsim_roi(mm,:,1) = mean(struct_m.res_xsim(:,2:n_reps,:),[1,2]);

    % Fill in the standard deviation for each one
    stat_rmse_roi(mm,:,2) = std(struct_m.res_rmse(:,2:n_reps,:),[],[1,2]);
    stat_xsim_roi(mm,:,2) = std(struct_m.res_xsim(:,2:n_reps,:),[],[1,2]);

    % % Fill in the mean and standard deviation for each one
    % gr_rmse_roi(mm,1) = mean(res_roi_rmse(:,2:n_reps),'all');
    % gr_rmse_roi(mm,2) = std(res_roi_rmse(:,2:n_reps),[],'all');
    % 
    % gr_xsim_roi(mm,1) = mean(res_roi_xsim(:,2:n_reps),'all');
    % gr_xsim_roi(mm,2) = std(res_roi_xsim(:,2:n_reps),[],'all');



end % for mm = 1:n_meth


%% Bar Chart of RMSE

% Figure
f1 = figure(1); clf;
set(f1,'Position',[150,300,(200 + 100.*n_rois),550]);

% Bar Chart
br1 = bar(stat_rmse_roi(:,:,1)',1,'FaceColor','Flat');
box on; hold on;

% Pre-allocate array of x positions for the ends of the boxes
xpos = zeros(n_rois,2);

for ee = 1:2
    xpos(:,ee) = br1(ee).XEndPoints;
end

% Pull out centre points and ends of the error bars
br_c = stat_rmse_roi(:,:,1)';
br_e = stat_rmse_roi(:,:,2)';

% Plot the error bars
er1 = errorbar(xpos(:),br_c(:),br_e(:),'k','LineStyle','none','LineWidth',1);

% Labels
ylim([0,135]);
ylabel('RMSE');
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
legend(names_methnice,'Location','NorthWest')
legend('boxoff');


%% Bar Chart of XSIM

% Figure
f2 = figure(2); clf;
set(f2,'Position',[250,300,(200 + 100.*n_rois),550]);

% Bar Chart
br2 = bar(stat_xsim_roi(:,:,1)',1,'FaceColor','Flat');
box on; hold on;

% Pre-allocate array of x positions for the ends of the boxes
xpos = zeros(n_rois,2);

for ee = 1:2
    xpos(:,ee) = br2(ee).XEndPoints;
end

% Pull out centre points and ends of the error bars
br_c = stat_xsim_roi(:,:,1)';
br_e = stat_xsim_roi(:,:,2)';

% Plot the error bars
er1 = errorbar(xpos(:),br_c(:),br_e(:),'k','LineStyle','none','LineWidth',1);

% Labels
ylim([0,1]);
ylabel('XSIM');
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
legend(names_methnice,'Location','NorthWest')
legend('boxoff');



%% ANOVA - Whole Brain

% % Choose data
% stat_table = stat_xsim_wbr(:,:)';
% 
% % Do one-way anova
% res_anova = anova(stat_table);
% 
% % Pull p-value
% res_pval = res_anova.stats.pValue(1);
% 
% % If there is a significance, do multiple comparisons
% if res_pval < 0.05
% 
%     t_mult_w = multcompare(res_anova);
% 
%     % Pull out only the statistically significant pairs
%     t_mult_w = t_mult_w(t_mult_w.pValue < 0.05, :);
% 
%     % Loop over remaining pairs
%     for pp = 1:size(t_mult_w,1)
%         disp([names_methnice{t_mult_w.Group1(pp)},' - ',names_methnice{t_mult_w.Group2(pp)},...
%              ' (WBR) are significantly different, at p = ',num2str(t_mult_w.pValue(pp))]);
%     end
% 
% 
% else % if res_pval < 0.05
%     disp('No significant differences (WBR)');
% 
% end % if res_pval < 0.05 // else
% 
% 
% %% ANOVA - ROIs
% 
% % Choose data
% stat_table = stat_xsim_roi(:,:)';
% 
% % Do one-way anova
% res_anova = anova(stat_table);
% 
% % Pull p-value
% res_pval = res_anova.stats.pValue(1);
% 
% % If there is a significance, do multiple comparisons
% if res_pval < 0.05
% 
%     t_mult_r = multcompare(res_anova);
% 
%     % Pull out only the statistically significant pairs
%     t_mult_r = t_mult_r(t_mult_r.pValue < 0.05, :);
% 
%     % Loop over remaining pairs
%     for pp = 1:size(t_mult_r,1)
%         disp([names_methnice{t_mult_r.Group1(pp)},' - ',names_methnice{t_mult_r.Group2(pp)},...
%              ' (ROI) are significantly different, at p = ',num2str(t_mult_r.pValue(pp))]);
%     end
% 
% 
% else % if res_pval < 0.05
%     disp('No significant differences (ROI)');
% 
% end % if res_pval < 0.05 // else


