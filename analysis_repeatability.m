% ANALYSIS_REPEATABILITY.m calculates repeatability statistics of QSM data.
%   Created by MT Cherukara, September 2024
%
%
%       Copyright (C) University College London, 2024
%
%
% QSM data should be reconstructed using QSM_PIPELINE.m and then ROI-average
% values extracted using ROI_AVERAGES.m. This script loads .mat files created by
% ROI_AVERAGES.m for the list of files specified by 'meth_names' and calculates
% ROI-based repeatability statistics for them.
%
% CHANGELOG:
%
% 2024-11-26 (MTC). Current version, with documentation.

clearvars;
close all;


%% Choose the Data

% Data directory
if isunix == 1
    dir_data = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/HNRepeatability_Data/';
else
    dir_data = '..\MTC_QSM\Analysis\HNRepeatability_Data\';
end

% Choose ROIs
% pickrois = [7,8,3,13:15];   % Neck
% pickrois = [9:12];            % Brain
pickrois = [9:12,7,8,3,13:15];  % All HN ROIs
n_pick = length(pickrois);

% Choose the data sets we want to load
% meth_names = {'unwrapped-SEGUE_bfr-PDF_susc-iterTik',...
%               'unwrapped-SEGUE_bfr-PDF_susc-StarQSM',...
%               'unwrapped-SEGUE_bfr-PDF_susc-FANSI',...
%               'unwrapped-SEGUE_bfr-PDF_susc-autoNDI',...
%               'unwrapped-SEGUE_bfr-PDF_susc-QSMnet',...
%               'unwrapped-SEGUE_bfr-TFI_susc-TFI'};
% meth_names = {'unwrapped-SEGUE_mask-nfv_bfr-LBV_susc-iterTik',...
%               'unwrapped-SEGUE_mask-nfv_bfr-PDF_susc-iterTik',...
%               'unwrapped-SEGUE_mask-nfv_bfr-VSHARP_susc-iterTik'};
meth_names = {'unwrapped-SEGUE_mask-ne_bfr-PDF_susc-autoNDI',...
              'unwrapped-SEGUE_mask-nfe_bfr-PDF_susc-autoNDI'};
n_meth = length(meth_names);

n_subs = 10;
n_reps = 6;

% Nice Method Names
% names_methnice = {'iterTik','StarQSM','FANSI','NDI','QSMnet','TFI'};
names_methnice = {'Noise Mask','Filled Noise Mask'};


% Bootstrap number of samples
n_boot = 1000;

% Do we want to plot bootstrap results
is_plotboot = 1;

% ROI degrees of freedom
roi_dof = [21,38,77,74,77,74,328,358,7490,2074,3222,1357,210,380,300]-1;

% For every statistic, we will store the actual value, followed by the lower and
% upper 95% confidence intervals
stat_rc = zeros(n_meth,n_pick,3);
stat_cv = zeros(n_meth,n_pick,3);
stat_ic = zeros(n_meth,n_pick,3);
stat_bsd = zeros(n_meth,n_pick,3);
stat_wsd = zeros(n_meth,n_pick,3);
stat_osd = zeros(n_meth,n_pick);
stat_mean = zeros(n_meth,n_pick);
stat_rsd = zeros(n_meth,n_pick,3);


%% Loop Through Methods and ROIs and Calculate Some Statistics

% Define some anonymous functions
std_intersub = @(x)(mean(std(x,[],1,'omitnan'),2,'omitnan'));
std_interrep = @(x)(mean(std(x,[],2,'omitnan'),1,'omitnan'));

% Pre-allocate main storage array
arr_av_full = zeros(n_subs.*n_reps,n_pick,n_meth);

for mm = 1:n_meth

    % Load the Data
    load(strcat(dir_data,'Repeatability_',meth_names{mm},'_ref_data.mat'));

    for rr = 1:length(pickrois)
    
        % Pull out this particular ROI
        r1 = pickrois(rr);
        arr_data = squeeze(arr_roi_av(:,:,r1));

        % Pull out the tissue standard deviation
        stat_rsd(mm,rr,1) = mean(arr_roi_sd(:,:,r1),'all','omitnan');

        % Calculate CIs on tissue standard deviation
        stat_rsd(mm,rr,2) = sqrt(roi_dof(r1)).*stat_rsd(mm,rr,1)./chi2inv(0.975,roi_dof(r1));
        stat_rsd(mm,rr,3) = sqrt(roi_dof(r1)).*stat_rsd(mm,rr,1)./chi2inv(0.025,roi_dof(r1));

        % Calculate within-subject standard deviation (and CIs)
        stat_wsd(mm,rr,1) = std_interrep(arr_data);
        stat_wsd(mm,rr,2:3) = bootci(n_boot,std_interrep,arr_data);
        % stat_wsd(mm,rr,4) = sqrt(roi_dof(r1)).*stat_wsd(mm,rr,1)./chi2inv(0.975,roi_dof(r1));
        % stat_wsd(mm,rr,5) = sqrt(roi_dof(r1)).*stat_wsd(mm,rr,1)./chi2inv(0.025,roi_dof(r1));

        % Calculate between-subject standard deviation (and CIs)
        stat_bsd(mm,rr,1) = std_intersub(arr_data);
        stat_bsd(mm,rr,2:3) = bootci(n_boot,std_intersub,arr_data);
        % stat_bsd(mm,rr,4) = sqrt(roi_dof(r1)).*stat_bsd(mm,rr,1)./chi2inv(0.975,roi_dof(r1));
        % stat_bsd(mm,rr,5) = sqrt(roi_dof(r1)).*stat_bsd(mm,rr,1)./chi2inv(0.025,roi_dof(r1));


        % Calculate overall mean
        stat_mean(mm,rr) = mean(arr_data,[1,2],'omitnan');

        % Calculate overall bias-corrected standard deviation
        stat_osd(mm,rr) = std(arr_data,[],[1,2],'omitnan').*n_reps./(n_reps-1);
   
        % Calculate Repeatability Coefficient (and Bootstrap)
        stat_rc(mm,rr,1) = mtc_repcoef(arr_data);
        stat_rc(mm,rr,2:3) = bootci(n_boot,@(x)[mtc_repcoef(x)],arr_data);
    
        % Calculate coefficient of variation
        stat_cv(mm,rr,1) = mtc_cofv(arr_data);
        stat_cv(mm,rr,2:3) = bootci(n_boot,@(x)[mtc_cofv(x)],arr_data);
    
        % Calculate ICC
        stat_ic(mm,rr,1) = mtc_icc(arr_data);
        stat_ic(mm,rr,2:3) = bootci(n_boot,@(x)[mtc_icc(x)],arr_data);

        % Store full data in a big array
        arr_av_full(:,rr,mm) = arr_data(:);
    
    end % for rr = 1:length(pickrois)

end % for mm = 1:n_meth



%% Box Plot of Actual Susceptibility Values

% % Make "grouping" variables
% data_group_rois = repmat(1:n_pick,[n_subs.*n_reps,1,n_meth]);
% data_group_meth = shiftdim(repmat((1:n_meth)',[1,n_subs.*n_reps,n_pick]),1);
% 
% % Create figure window
% f8 = figure(8); clf;
% set(f8,'Position',[250,300,(200 + 20.*n_meth.*n_pick),550]);
% 
% % Box chart
% bx8 = boxchart(1.25.*data_group_rois(:),arr_av_full(:),'GroupByColor',data_group_meth(:),...
%                'BoxWidth',0.8);
% 
% % Straight line at 0
% hold on;
% plot([0.5,(1.25.*n_pick)+0.75],[0,0],'k:','LineWidth',1.5);
% 
% % Labels and decoration
% box on;
% xticks(1.25.*(1:n_pick));
% xticklabels(names_roi(pickrois));
% xlim([0.5,(1.25.*n_pick)+0.75]);
% legend(names_methnice,'Location','NorthWest');
% legend('boxoff');
% ylabel('Susceptibility (ppm)');
% set(gca,'FontSize',16,'FontName','Calibri');


%% Compare Standard Deviations 

% Pull out the three standard deviations
mat_sd = [stat_rsd(:,:,1); stat_wsd(:,:,1); stat_bsd(:,:,1)]';
mat_sdlb = [stat_rsd(:,:,2); stat_wsd(:,:,2); stat_bsd(:,:,2)]';
mat_sdub = [stat_rsd(:,:,3); stat_wsd(:,:,3); stat_bsd(:,:,3)]';

% x-axis positions
xpos = repmat(-0.25:0.1:0.25,n_pick,1) + repmat((1:n_pick)',1,6);

% Create figure window
f9 = figure(11); clf;
set(f9,'Position',[250,350,(200 + 100.*n_pick),550]);

% Fake extra error bar plots which won't be seen but are just used to create
% extra terms for the legend
er1 = plot([1,1],[5,5],'LineStyle','None','Marker','x','LineWidth',1.5,'MarkerSize',10);
hold on; box on;
er2 = plot([2,2],[5,5],'LineStyle','None','Marker','o','LineWidth',1.5);
er3 = plot([3,3],[5,5],'diamond','LineStyle','None','Marker','diamond','LineWidth',1.5);
er1.Color = [0, 0, 0];
er2.Color = [0, 0, 0];
er3.Color = [0, 0, 0];

% Plot intervals as error bars
er9 = errorbar(xpos,mat_sd,mat_sdlb,mat_sdub,'Marker','x',...
              'LineStyle','none','LineWidth',1.5);

% Manually Colour and Texture the Bars
er9(1).Color = [0.000, 0.447, 0.741];
er9(2).Color = [0.850, 0.325, 0.098];
er9(3).Color = [0.000, 0.447, 0.741];
er9(4).Color = [0.850, 0.325, 0.098];
er9(5).Color = [0.000, 0.447, 0.741];
er9(6).Color = [0.850, 0.325, 0.098];

er9(3).Bar.LineStyle = 'dashed';
er9(4).Bar.LineStyle = 'dashed';
er9(5).Bar.LineStyle = 'dotted';
er9(6).Bar.LineStyle = 'dotted';

er9(1).MarkerSize = 10;
er9(2).MarkerSize = 10;
er9(3).Marker = 'o';
er9(4).Marker = 'o';
er9(5).Marker = 'diamond';
er9(6).Marker = 'diamond';

% Labels
ylabel('Standard Deviation (ppm)');
xticks(1:n_pick);
xticklabels(names_roi(pickrois));
legend([{'Within-Region SD'},{'Within-Subject SD'},{'Between-Subject SD'},names_methnice],'Location','NorthWest');
legend('boxoff');
set(gca,'FontSize',16,'FontName','Calibri');
ylim([0,0.22]);
xlim([0.25,10.75]);


%% Plot Standard Deviation

f10 = figure(10); clf;
set(f10,'Position',[100,100,(300 + 35.*n_meth.*n_pick),550]);

b10 = bar(stat_wsd(:,:,1)',1,'FaceColor','Flat');
box on; hold on;

% Plot error bars showing Bootstrap results
if is_plotboot == 1

    % Pre-allocate array of x positions for the ends of the boxes
    xpos = zeros(n_pick,n_meth);

    % Loop through methods and fill it in
    for mm = 1:n_meth
        xpos(:,mm) = b10(mm).XEndPoints;
    end

    % Pull out centre points and ends of the error bars
    eb_mid = stat_wsd(:,:,1)';
    eb_bot = eb_mid - stat_wsd(:,:,2)';
    eb_top = stat_wsd(:,:,3)' - eb_mid;

    % Plot the error bars
    er10 = errorbar(xpos(:),eb_mid(:),eb_bot(:),eb_top(:),...
                    'k','LineStyle','none','LineWidth',1);
end

% Labels
ylabel('Within-Subject SD (ppm)');
ylim([0,0.1])
xticks(1:n_pick);
xticklabels(names_roi(pickrois));
legend(names_methnice,'Location','NorthWest');
legend('boxoff');
set(gca,'FontSize',16,'FontName','Calibri');


%% Plot Repeatability Coefficient

% % Create Figure Window
% f11 = figure(21); clf;
% set(f11,'Position',[100,100,(300 + 35.*n_meth.*n_pick),550]);
% 
% % Bar Chart
% b11 = bar(stat_rc(:,:,1)',1,'FaceColor','Flat','FaceAlpha',1);
% box on; hold on;
% 
% % Plot error bars showing Bootstrap results
% if is_plotboot == 1
% 
%     % Pre-allocate array of x positions for the ends of the boxes
%     xpos = zeros(n_pick,n_meth);
% 
%     % Loop through methods and fill it in
%     for mm = 1:n_meth
%         xpos(:,mm) = b11(mm).XEndPoints;
%     end
% 
%     % Pull out centre points and ends of the error bars
%     eb_mid = stat_rc(:,:,1)';
%     eb_bot = eb_mid - stat_rc(:,:,2)';
%     eb_top = stat_rc(:,:,3)' - eb_mid;
% 
%     % Plot the error bars
%     er11 = errorbar(xpos(:),eb_mid(:),eb_bot(:),eb_top(:),...
%                     'k','LineStyle','none','LineWidth',1);
% end
% 
% % Labels
% ylabel('Repeatability Coefficient (ppm)');
% % ylim([0,0.25])
% xticks(1:n_pick);
% xticklabels(names_roi(pickrois));
% legend(names_methnice,'Location','NorthWest');
% legend('boxoff');
% set(gca,'FontSize',16,'FontName','Calibri');



%% Plot Coefficient of Variation

% % Create Figure Window
% f12 = figure(22); clf;
% set(f12,'Position',[100,100,(300 + 35.*n_meth.*n_pick),550]);
% 
% % Bar Chart
% b12 = bar(stat_cv(:,:,1)',1,'FaceColor','Flat','FaceAlpha',1);
% box on; hold on;
% 
% % Lines
% xv = [0.45,n_pick+0.55];
% plot(xv,[0.2,0.2],'k:','LineWidth',1.5);
% % plot(xv,[0.15,0.15],'k--');
% 
% % Labels
% ylabel('Coefficient of Variation');
% ylim([0,2])
% xlim([0.5,n_pick+0.5])
% xticks(1:n_pick);
% xticklabels(names_roi(pickrois));
% legend(names_methnice,'Location','NorthWest');
% % legend('boxoff');
% set(gca,'FontSize',16,'FontName','Calibri');


%% Plot ICC

% Create Figure Window
f13 = figure(13); clf;
set(f13,'Position',[200,100,(300 + 35.*n_meth.*n_pick),550]);

% Bar Chart
b13 = bar(stat_ic(:,:,1)',1,'FaceColor','Flat','FaceAlpha',1);
box on; hold on;

% Lines
xv = [0.45,n_pick+0.55];
plot(xv,[0.75,0.75],'k:','LineWidth',1.5);
plot(xv,[0.9,0.9],'k--','LineWidth',1.25);

% Plot error bars showing Bootstrap results
if is_plotboot == 1

    % Pre-allocate array of x positions for the ends of the boxes
    xpos = zeros(n_pick,n_meth);

    % Loop through methods and fill it in
    for mm = 1:n_meth
        xpos(:,mm) = b13(mm).XEndPoints;
    end

    % Pull out centre points and ends of the error bars
    eb_mid = stat_ic(:,:,1)';
    eb_bot = eb_mid - stat_ic(:,:,2)';
    eb_top = stat_ic(:,:,3)' - eb_mid;

    % Plot the error bars
    er11 = errorbar(xpos(:),eb_mid(:),eb_bot(:),eb_top(:),...
                    'k','LineStyle','none','LineWidth',1);
end

% Labels
ylabel('Intra-class Correlation Coefficient');
ylim([0,1])
xlim([0.5,n_pick+0.5])
xlim([0.5,n_pick+0.5])
xticks(1:n_pick);
xticklabels(names_roi(pickrois));
legend(names_methnice,'Location','SouthEast');
% legend('boxoff');
set(gca,'FontSize',16,'FontName','Calibri');


