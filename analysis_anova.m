% ANALYSIS_ANOVA.m performs repeated-measurs ANOVA on QSM data.
%   Created by MT Cherukar, September 2024
%
%
%       Copyright (C) University College London, 2024
%
% 
% Load some ROI-averaged susceptibility data (calculated by ROI_AVERAGES.m) and
% performs Repeated Measures Analysis of Variation (RANOVA) on it. Displays the
% results as a nice looking Heatmap, and stores them as a Table (you will need
% to save this manually to keep it).
%
% CHANGELOG:
%
% 2024-11-26 (MTC). Original version, with documentation.


clearvars;
close all;

% mtc_setup;


%% Choose the Data

% Data directory
if isunix == 1
    dir_data = '/home/cherukara/Documents/Coding/MTC_QSM/Analysis/HNRepeatability_Data/';
else
    dir_data = '..\MTC_QSM\Analysis\HNRepeatability_Data\';
end

% Choose ROIs
% pickrois = [3,7,8,13:15];   % Neck
pickrois = 9:12;            % Brain
n_pick = length(pickrois);

% Choose the data sets we want to load
meth_names = {'unwrapped-SEGUE_bfr-PDF_susc-iterTik',...
              'unwrapped-SEGUE_bfr-PDF_susc-StarQSM',...
              'unwrapped-SEGUE_bfr-PDF_susc-FANSI',...
              'unwrapped-SEGUE_bfr-PDF_susc-autoNDI',...
              'unwrapped-SEGUE_bfr-PDF_susc-QSMnet',...
              'unwrapped-SEGUE_bfr-TFI_susc-TFI'};
n_meth = length(meth_names);

% Nice Method Names
names_methnice = {'iterTik','StarQSM','FANSI','autoNDI','QSMnet','TFI'};

% Lengths
n_subs = 10;
n_reps = 6;


%% Load everything into a big array

% Pre-allocate the array
%       ARR_DATA_ALL has dimensions SUBS * REPEATS * ROIS * METHODS
arr_data_all = zeros(n_subs,n_reps,n_pick,n_meth);

% Loop over Methods
for mm = 1:n_meth

    % Load in the data
    %       ARR_ROI_AV has dimensions SUBS * REPEATS * ROIS
    load(strcat(dir_data,'Repeatability_',meth_names{mm},'_data.mat'));

    % Pick out the data we want and stick it in a big array
    arr_data_all(:,:,:,mm) = arr_roi_av(:,:,pickrois);

end % for mm = 1:n_meth


%% Pre-allocate storage arrays and define pairs and names

% Pre-allocate results arrays
%   First column is the actual statistic, second column is the p-value
res_mauch = zeros(n_pick,2);
res_ranova = zeros(n_pick,2);
res_single = zeros(n_pick,n_meth,2);

% This will be for storing the pairwise MULTCOMPARE results,
%   First lower bound, then difference, then upper bound, then p-value
res_mult = zeros(n_pick,nchoosek(n_meth,2),4);

% Pre-allocate a cell array to put all the tables in
arr_roitables = cell(1,n_pick);

% List the actual pairs
pairs_mult = nchoosek(1:n_meth,2);

% Add a fake one in to PAIRS_MULT
pairs_mult = [pairs_mult;[n_meth,1]];
n_pair = size(pairs_mult,1);

% Pre-allocate cell array
names_pairs = cell(1,n_pair);

% Make pair names
for pp = 1:n_pair

    names_pairs{pp} = [names_methnice{pairs_mult(pp,1)},'-',names_methnice{pairs_mult(pp,2)}];

end % for pp = 1:n_pair

% Create an array of subject IDs (this will be a dummy variable in the ANOVA)
str_subs = arrayfun(@(x) sprintf('s%d', x), 1:n_subs, 'UniformOutput', false)';
arr_subs = repmat(str_subs,1,n_reps);


%% Do ANOVA analysis for each ROI


% Loop over ROIs
for rr = 1:n_pick

    % Pull out an array for this ROI
    %       ARR_DATA_ROI has dimensions SUBS * REPEATS * METHODS
    arr_data_roi = squeeze(arr_data_all(:,:,rr,:));

    % Reshape
    %       ARR_TABLE_ROI has dimensions (SUBS * REPEATS) * METHODS
    arr_table_roi = reshape(arr_data_roi,[n_subs.*n_reps,n_meth]);

    % Make it a table (this is currently hardcoded to 6 methods)
    t_roidata = table(arr_subs(:),arr_table_roi(:,1),arr_table_roi(:,2),...
                                  arr_table_roi(:,3),arr_table_roi(:,4),...
                                  arr_table_roi(:,5),arr_table_roi(:,6),...
                      'VariableNames',{'scan','meas1','meas2','meas3','meas4','meas5','meas6'});

    % Make a single-row WithinDesign table
    t_single = table(1,'VariableNames',{'Measurements'});

    % Loop over Methods and do single-method ANOVA for each one
    for mm = 1:n_meth

        % Design string
        str_design = ['meas',num2str(mm),'~scan'];

        % Fit a repeated measures model within a single method ("Measurement")
        rm_single = fitrm(t_roidata,str_design,'WithinDesign',t_single);
    
        % Do ANOVA of this data (not RANOVA, because we only have one measure)
        t_ranova1 = anova(rm_single);
    
        % Store the ANOVA results
        %   The Constant:scan term tells us if there is a significant difference
        %   between subjects ("scans") for this method ("Measurement"), so we want
        %   to store that F statistic and pValue
        res_single(rr,mm,1) = t_ranova1.F(2);
        res_single(rr,mm,2) = t_ranova1.pValue(2);

    end % for mm = 1:n_meth

    % Make a WithinDesign table for all the different methods
    t_meas = table((1:n_meth)','VariableNames',{'Measurements'});

    % Fit a repeated measures model across all measurements
    rm_all = fitrm(t_roidata,'meas1-meas6~scan','WithinDesign',t_meas);

    % Do Mauchly test
    t_mauch = mauchly(rm_all);

    % Store Mauchly test results
    %   We need to store the Chi-statistic and the p-value
    res_mauch(rr,:) = [t_mauch.ChiStat(1), t_mauch.pValue(1)];

    % Do repeated measures ANOVA
    t_ranova = ranova(rm_all);

    % Store ANOVA results
    %   The most important thing is the (Intercept):Measurements, since this
    %   tells us if there is a significant difference between Measurements
    %   (which are the methods we are comparing).
    %   
    %   The other field is scan:Measurements, which tells us if there is a
    %   significant variation across subjects and different methods. This is not
    %   necessarily interesting.
    %
    %   We want to store the F-statistic and the p-value of the (Intercept):
    %   Measruements effect 
    res_ranova(rr,1) = t_ranova.F(2);

    if res_mauch(rr,2) < 0.05
        % If Mauchly's test is violated, use the lower bound on the p-value
        res_ranova(rr,2) = t_ranova.pValueLB(2);
    else
        % Otherwise, use the normal ANOVA p-value
        res_ranova(rr,2) = t_ranova.pValue(2);
    end
    
    % Make multiple comparisons between Methods
    t_mult = multcompare(rm_all,'Measurements');

    % Make an empty cell array
    arr_roidiffs = cell(n_pair,6);

    % Store the pairwise results

    % Loop over comparison pairs
    for pp = 1:n_pair

        % Find the correct row corresponding to the current pair
        i_row = find(t_mult.Measurements_1 == pairs_mult(pp,1) & ...
                     t_mult.Measurements_2 == pairs_mult(pp,2) );

        % Store the names
        arr_roidiffs{pp,1} = string(names_methnice(pairs_mult(pp,1)));
        arr_roidiffs{pp,2} = string(names_methnice(pairs_mult(pp,2)));

        % Store the difference and its upper and lower bounds (converted to ppb)
        arr_roidiffs{pp,3} = 1000.*abs(t_mult.Difference(i_row));
        arr_roidiffs{pp,4} = 1000.*min(abs(t_mult.Lower(i_row)),abs(t_mult.Upper(i_row)));
        arr_roidiffs{pp,5} = 1000.*max(abs(t_mult.Lower(i_row)),abs(t_mult.Upper(i_row)));

        % Store the p-value of that row
        arr_roidiffs{pp,6} = t_mult.pValue(i_row);

        % Store the difference, and its upper and lower bounds
        res_mult(rr,pp,1) = t_mult.Lower(i_row);
        res_mult(rr,pp,2) = t_mult.Difference(i_row);
        res_mult(rr,pp,3) = t_mult.Upper(i_row);

        % Store the p-value of that row
        res_mult(rr,pp,4) = t_mult.pValue(i_row);
        
    end % for pp = 1:n_pair

    % % Make a table for storing the results
    % t_roi = table(1000.*res_mult(rr,:,2)',1000.*[res_mult(rr,:,1)',res_mult(rr,:,3)'],...
    %               res_mult(rr,:,4)' ,res_mult(rr,:,4)' < 0.05,...
    %               'VariableNames',{'Difference','95% CI','pValue','Significant?'},'RowNames',names_pairs);

    t_roi = cell2table(arr_roidiffs,'VariableNames',{'Method 1','Method 2','Difference','Lower','Upper','pValue'});

    % Store the table in a cell array
    arr_roitables{rr} = t_roi;

end % for rr = 1:n_pick


%% Store Results in a nice Table

% Pre-allocate cell array
names_pairs = cell(1,n_pair);

% Make pair names
for pp = 1:n_pair

    names_pairs{pp} = [names_methnice{pairs_mult(pp,1)},'-',names_methnice{pairs_mult(pp,2)}];

end % for pp = 1:n_pair

% Create table of multiple-comparisons results
t_multres = array2table(squeeze(res_mult(:,:,4) < 0.05),'VariableNames',names_pairs,'RowNames',names_roi(pickrois));

% Make a table of the single-method ANOVA results
t_single = array2table(squeeze(res_single(:,:,2) < 0.05),'VariableNames',names_methnice,'RowNames',names_roi(pickrois));


%% Display Heatmap

% Loop over ROIs
for rr = 1:n_pick

    % Pull out the current table
    t_roi = arr_roitables{rr};

    % Fake the last variable
    t_roi.Difference(end) = NaN;

    % Make a figure
    f00 = figure(21+rr); clf;
    set(f00,'Position',[300,200+(50.*rr),600,500]);

    % Plot heatmap
    h00 = heatmap(t_roi,'Method 2','Method 1','ColorVariable','Difference');

    % Labels
    h00.XDisplayData = names_methnice;
    h00.YDisplayData = names_methnice;
    h00.XLabel = '';
    h00.YLabel = '';
    h00.MissingDataLabel = '';
    h00.ColorLimits = [0,200];
    h00.Title = names_roi(pickrois(rr));
    set(gca,'FontSize',16,'FontName','Calibri');
    set(gca,'ColorbarVisible','Off');

    % Save Figure
    exportgraphics(f00,sprintf('Heatmap_RANOVA_%d_%s.tif',rr,names_roi{pickrois(rr)}));


end % for rr = 1:n_pick