%% CAPSTEP: TEP-MEP RELATIONSHIP
% Written by Dominika for the CAPS-TEP project (2022)

% ----- Prepare data - matched single trials ----- 

%% parameters
clear all
clc

% choose the Git folder --> source of saved default files
path = 'C:\Users\sulcova\Desktop\GitHub\CAPS-TEP';
folder_git = uigetdir(path, 'Choose the Git folder');

% choose the folder with processed data 
path = 'E:\Data\CAPS-TEP - data\Processed data';
folder_input = uigetdir(path, 'Coose the input folder');

% choose the results and figures folder  
path = 'E:\UCL\O365G-NOCIONS - CAPS-TEP - CAPS-TEP\Results';
folder_results = uigetdir(path, 'Coose the results folder');
folder_figures = [folder_results '\CAPSTEP_figures'];

% dataset
participant = [1:10, 12:14, 16:18, 20, 22:24];
condition = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix_TEP = 'chanflip avg bl icfilt ica visual crop notch bandpass prefilt prea ds art-sup dc ep reref chan-select reref EEG';

% visualization
load([folder_git '\colours.mat'])
figure_counter = 1;

% output 
output_file = [folder_results '\CAPSTEP_output.mat'];

%% TEP-MEP CORRELATION - individual average
% ----- section input -----
outliers = [9, 13, 23];
remove_outliers = 0;
TEP_peak = 6;
% -------------------------
% load data
load(output_file, 'CAPSTEP_TEP_peaks', 'CAPSTEP_MEP', 'CAPSTEP_TEP_default') 

% organize data
for c = 1:length(condition)
    for t = 1:length(time)
        for p = 1:length(participant)
            % determine the subject
            if participant(p) < 10
               subject = ['0' num2str(participant(p))];
            else
               subject = num2str(participant(p)); 
            end
            
            % determine row
            row = (c-1)*length(time)*length(participant) + (t-1)*length(participant) + p;
            
            % TEP peak amplitude = x 
            x(row) = CAPSTEP_TEP_peaks.amplitude_peak(c, t, p, TEP_peak);

            % MEP p2p amplitude = y
            y(row) = CAPSTEP_MEP.amplitude_zero(string(CAPSTEP_MEP.subject) == subject &...
                string(CAPSTEP_MEP.session) == condition{c} &...
                string(CAPSTEP_MEP.timepoint) == time{t}); 
            
            % choose point colours
            if c == 1
                marker_col(row, :) = [1,0.68,0.68];
            elseif c == 2
                marker_col(row, :) = [0.69,0.87,1];
            end
        end
    end
end
clear c t p row

% rank the data
[temp, x_ranked]  = ismember(x, unique(x));
[temp, y_ranked]  = ismember(y, unique(y));
clear a temp

% prepare linear model: y ~ 1 + x
data_model = fitlm(x_ranked, y_ranked, 'VarNames', {['TEP_' CAPSTEP_TEP_default.peaks{TEP_peak}] 'MEP'});
    
% extract statistics
data_table = data_model.Coefficients;
data_table.R2(1) = data_model.Rsquared.Ordinary;
data_table.R2(2) = data_model.Rsquared.Adjusted;
data_table.Properties.RowNames = {};
data_table.Variable(1) = {'(intercept)'};
data_table.Variable(2) = data_model.VariableNames(1);

% plot data + regression line
if CAPSTEP_TEP_default.peaks{TEP_peak}(1) == 'N'
    flip_x = true;
else
    flip_x = false;
end
fig = figure(figure_counter);
hold on
plot_corr(data_model, cat(1, x_ranked, y_ranked)', marker_col, 'Spearman', flip_x)

% name and save figure
figure_name = sprintf('CAPSTEP_TEPxMEP_average_%s', CAPSTEP_TEP_default.peaks{TEP_peak});
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])

% update figure counter
figure_counter = figure_counter + 1 ;

%% functions
function plot_corr(data_model, data_corr, marker_col, corr_type, flip_x)
% plot correlation
plot_cor = plotAdded(data_model);

% identify p    
[cor_coef, cor_p] = corr(data_corr, 'Type', corr_type);

% adjust parameters    
set(gca, 'FontSize', 20)
xlabel(data_model.VariableNames{1}); 
ylabel(data_model.VariableNames{2});
plot_cor(2).Color = [0 0 0]; plot_cor(2).LineWidth = 4; 
plot_cor(3).Color = [0 0 0]; plot_cor(3).LineWidth = 2;
legend off
title('')

% replot markers
for c = 1:size(data_corr, 1)
    scatter(data_corr(c, 1), data_corr(c, 2), 90, marker_col(c, :), 'filled');
    hold on
end

% flip x if required
if flip_x
    set(gca, 'XDir','reverse')
end

% adjust ticks
xticklabels({'most negative','', '', 'most positive'})
yticklabels({'low','', '', 'high'})

% add annotations
text_pos = [0.90 0.78 0.66];
% rectangle('Position', [2, 27, 18, 12], 'FaceColor',[1 1 1], 'EdgeColor', [0.5 0.5 0.5])
T(1) = text(0.1, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.1, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Adjusted), 'Units', 'Normalized');
T(3) = text(0.1, text_pos(3), sprintf('p = %1.5f', cor_p(1, 2)), 'Units', 'Normalized');
set(T(1), 'fontsize', 20, 'fontangle', 'italic'); 
set(T(2), 'fontsize', 20, 'fontweight', 'bold'); 
set(T(3), 'fontsize', 20, 'fontweight', 'bold'); 

end