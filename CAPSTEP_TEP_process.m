%% CAPSTEP: PROCESS TEPs
% Written by Dominika for the CAPS-TEP project (2022)

% ----- Extract peak values of TEP components ----- 
% 1) load data
%     - data cropped [-0.05 0.3]s 
%     --> datasets saved to structure 'CAPSTEP_TEP_data' and added
%     to the global MATLAB output file
% 2) calculate GFP, extract peak latencies
%     - calculate average across all datasets and corresponding GFP
%     - plot GFP and mark peaks ('TEP_GFP')
%     --> mean data saved in 'CAPSTEP_TEP_mean' and added to the output file
% 3) identify EOIs
%     - based on class map voltages obtained by MS analysis
%     - frontal and temporal electrodes excluded due to possible noisiness
%     --> EOI labels saved to structure 'CAPSTEP_TEP_default' and added
%     to the global MATLAB output file
% 4) extract TEP peak amplitudes
%     --> peak values saved to structure 'CAPSTEP_results' and added tp the
%     global MATLAB output file
% 5) export for R
%     - save in a long-format table suitable for statistics in R
%     --> 'CAPSTEP_TEP_values.csv'

% ----- Group visualization ----- 
% 6) plot baseline TEPs
% 7) plot peak values - boxplot
% 8) plot ratings - lineplot
% 9) plot peak values + ratings - lineplot

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
prefix = 'chanflip avg bl icfilt ica visual crop notch bandpass prefilt prea ds art-sup dc ep reref chan-select reref EEG';
condition = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};

% output 
output_file = [folder_results '\CAPSTEP_output.mat' ];

% % TEP default
% CAPSTEP_TEP_default.electrodes = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2','CP1','CP2','FC5','FC6','CP5','CP6','C1','C2'};
% CAPSTEP_TEP_default.peaks = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'}; 
% save(output_file, 'CAPSTEP_TEP_default', '-append');

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.15;

% visualization 
figure_counter = 1;
load([folder_input '\' prefix ' 10 caps baseline.lw6'], '-mat')
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for c = 1:length(condition) 
            for t = 1:length(time)
               colours(c, (t-1)*3 + 1:(t-1)*3 + 3) = uisetcolor; 
            end
        end
        save([folder_git '\colours.mat'], 'colours'); 
    case 'NO'
    load([folder_git '\colours.mat'])
end
clear answer c t 

%% 1) extraction of individual data
% load data
for c = 1:length(condition)
    for t = 1:length(time)
        for p = 1:length(participant)
            % define participant
            if participant(p) < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end
            
            % load dataset
            load([folder_input '\' prefix ' ' subj ' ' condition{c} ' ' time{t} '.mat'])

            % append the data in the data matrix
            CAPSTEP_TEP_data(c, t, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                  
        end
    end
end
clear c t p subj data
disp(['Data import finished. Datasize: ' num2str(size(CAPSTEP_TEP_data))])

% save to the global MATLAB file
save(output_file, 'CAPSTEP_TEP_data', '-append');

%% 2) GFP - peak identification
% ----- section input -----
labeled = 'off';
max_peaks = 6;
% participant = [];
% -------------------------
% load data if necessary
if exist('CAPSTEP_TEP_data') ~= 1
    load(output_file, 'CAPSTEP_TEP_data', 'CAPSTEP_TEP_default')
end

% calculate mean values 
CAPSTEP_TEP_mean.data = squeeze(mean(CAPSTEP_TEP_data, [1 2 3]));
CAPSTEP_TEP_mean.GFP = std(CAPSTEP_TEP_mean.data, 1);

% launch the figure
fig = figure(figure_counter);
hold on

% extract peak latencies
h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
TEP_peaks = gfp_plot(x, CAPSTEP_TEP_mean.GFP, time_window, xstep, labeled, 'max_peaks', max_peaks);
title('GFP of grand average TEPs', 'fontsize', 16, 'fontweight', 'bold')

% encode peak center --> add N45 if not identified
CAPSTEP_TEP_default.center = [TEP_peaks(1:2), 0.045, TEP_peaks(3:end)];

% choose data for topoplots 
for e = 1:size(CAPSTEP_TEP_mean.data, 1)
    for i = 1:size(CAPSTEP_TEP_mean.data, 2)
        data_topoplot(1, e, 1, 1, 1, i) = CAPSTEP_TEP_mean.data(e, i);
    end
end
clear e i

% add topoplots
for k = 1:length(TEP_peaks)
    % plot the topoplot
    h_axis(1 + k) = subplot(3, max_peaks, 2*max_peaks + k);
    topo_plot(header, data_topoplot, TEP_peaks(k), time_window(1), [-3, 3])

    % shift down
    pos = get(h_axis(1 + k), 'Position');
    pos(2) = pos(2) - 0.05;
    set(h_axis(1 + k), 'Position', pos);

    % add timing
    text(-0.3, -0.8, sprintf('%1.0f ms', TEP_peaks(k)*1000), 'Color', [1 0 0], 'FontSize', 14)
end
hold off

% save figure
savefig([folder_figures '\CAPSTEP_TEP_GFP.fig'])
saveas(fig, [folder_figures '\CAPSTEP_TEP_GFP.png'])

% update figure counteer
figure_counter = figure_counter + 1;

% append new variable to the general MATLAB file
save(output_file, 'CAPSTEP_TEP_mean', 'CAPSTEP_TEP_default', '-append');

clear c t k data_topoplot fig figure_name pos h_axis TEP_peaks labeled max_peaks 

%% 3) identify EOIs
% ----- section input -----
eoi_n = [3,1,3,3,3,3];
% participant = [];
% -------------------------
% load data if necessary
if exist('CAPSTEP_TEP_data') ~= 1
    load(output_file, 'CAPSTEP_TEP_data', 'CAPSTEP_TEP_default')
end

% extract voltage from MS maps
for k = 1:length(CAPSTEP_TEP_default.peaks)
    MS_voltage(k, :) = readmatrix([folder_results '\CAPSTEP_TEP_analysis\MS_' num2str(k) '.txt']); 
end 
clear k

% choose EOIs based on the voltage
for k = 1:length(CAPSTEP_TEP_default.peaks)
    % index max/min values
    if strcmp(CAPSTEP_TEP_default.peaks{k}(1), 'P')
        [eoi_val, eoi_i] = maxk(MS_voltage(k, :), eoi_n(k));
    else
        [eoi_val, eoi_i] = mink(MS_voltage(k, :), eoi_n(k));
    end
    
    % identify EOI labels
    CAPSTEP_TEP_default.eoi{k} = CAPSTEP_TEP_default.electrodes(eoi_i);
    disp([CAPSTEP_TEP_default.peaks{k} ': ' CAPSTEP_TEP_default.eoi{k}])
end
clear k eoi_val eoi_i eoi_n

% save to the global MATLAB file
save(output_file, 'CAPSTEP_TEP_default', '-append');

%% 4) TEP amplitude extraction
% ----- section input -----
CAPSTEP_TEP_default.span = [0.014 0.02 0.014 0.04 0.05 0.06];
percent = 20;                           % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                      % y limits for topoplots 
% participant = [];
% -------------------------
% set colours for visualisation
col_fig = [1.0000    0.4118    0.1608; 0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.7176    0.2745    1.0000; 0.3647    0.2078    0.9882];

% loop through subjects
for p = 1:length(participant)
    % setup names   
    figure_title = sprintf('Subject n. %d', participant(p));    
       
    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;
        
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
    % loop through peaks
    for k = 1:length(CAPSTEP_TEP_default.peaks)   
        % identify peak polarity
        if strcmp(CAPSTEP_TEP_default.peaks{k}(1), 'P')
            polarity = 1;
        else
            polarity = -1;
        end
            
        % identify EOIs
        eoi = [];
        for e = 1:length(CAPSTEP_TEP_default.eoi{k})
            eoi(e) = find(strcmp(CAPSTEP_TEP_default.electrodes, CAPSTEP_TEP_default.eoi{k}{e}));
        end
            
        % choose data
        for c = 1:length(condition)
            for t = 1:length(time)
                % data for timeseries visualization
                data_visual((c-1)*length(time) + t, :, :) = squeeze(mean(CAPSTEP_TEP_data(c, t, p, eoi, :), 4)); 
                
                % data for topoplot
                for e = 1:30
                    data_topoplot((c-1)*length(time) + t, e, 1, 1, 1, :) = squeeze(CAPSTEP_TEP_data(c, t, p, e, :));
                end
            end
        end
                       
        % define default TOI 
        center = CAPSTEP_TEP_default.center(k);
        span = CAPSTEP_TEP_default.span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;                
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 6, 1:18);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, CAPSTEP_TEP_default.peaks{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 6, 1:18)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % extract amplitude and latency
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, time_window(1), polarity);

            % update the figure            
            subplot(4, 6, 1:18)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
            end
            plot(x, data_visual(1:size(data_visual, 1)/2, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 2.5)
            plot(x, data_visual(size(data_visual, 1)/2 + 1:end, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5)
            plot(lat_peak, y_max, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')

            % add the topoplot   
            b = 1;
            for a = [1, size(data_visual, 1)/2, size(data_visual, 1)/2+1, size(data_visual, 1)]
                subplot(4, 6, 18 + b);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims) 
                b = b + 1;
            end
            clear a b

            % ask for approval
            answer = questdlg('Do you want to proceed?', CAPSTEP_TEP_default.peaks{k},...
                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');

            % switch action
            switch answer
                case 'Yes, extract outcome values.'
                    % close the figure
                    close(fig_1)

                    % exit the while loop
                    finish = 1;

                case 'No, I will adjust TOI manually.'
                    % assign previous center and span
                    choose_center = center;  
                    choose_span = 2 * span;  

                    % identify the limits for visualisation of current peak
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / xstep);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / xstep);
                    choose_x = (choose_center - choose_span/2) : xstep : (choose_center + choose_span/2);

                    % prepare data and header for visualization
                    choose_data = data_visual(:, choose_x1 : choose_x2);
                    choose_header = header;
                    choose_header.datasize(6) = length(choose_data);  
                    choose_header.xstart = choose_center - choose_span/2;

                    % check if vector size matches
                    if size(choose_data, 2) ~= length(choose_x)
                        diff = size(choose_data, 2) - length(choose_x);
                        if diff > 0
                            choose_data = choose_data(:, 1:end - diff);
                        elseif diff < 0
                            choose_x = choose_x(1:end + diff);
                        end
                    end

                    % launch the choosing figure                 
                    choose_figure_name = ['Choose manually peak ' CAPSTEP_TEP_default.peaks{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2, 'Color', [0.0549    0.5216    0.8118])
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    for a = 1
                        choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims);
                    end           

                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    for a = 1
                        subplot(4, 4, a) 
                        cla(choose_axesHandles(2))
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims);
                    end
                    hold off

                    % update the central latency
                    center = pos_x;

                    % close the choosing figure
                    pause(2)
                    close(choose_fig)

                    % close the the main figure
                    close(fig_1)
                end
        end
        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
            choose data choose_center choose_axesHandles answer

        % record outcome variables
        for c = 1:length(condition)
            for t = 1:length(time)
                CAPSTEP_TEP_peaks.latency(c, t, p, k) = lat_peak((c-1)*length(time) + t); 
                CAPSTEP_TEP_peaks.amplitude_peak(c, t, p, k) = y_max((c-1)*length(time) + t); 
                CAPSTEP_TEP_peaks.amplitude_mean(c, t, p, k) = y_mean((c-1)*length(time) + t); 
            end
        end

        % set up the main figure
        figure(fig)
        subplot(length(CAPSTEP_TEP_default.peaks), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', CAPSTEP_TEP_default.peaks{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = {'CAPS - baseline' 'CAPS - last' 'CTRL - baseline' 'CTRL - last'};
        b = 1;
        for a = [1, size(data_visual, 1)/2, size(data_visual, 1)/2+1, size(data_visual, 1)]
            subplot(length(CAPSTEP_TEP_default.peaks), 7, 3 + b + 7*(axis_counter-1))
            topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims);   
            if axis_counter == 1
                 set(get(gca, 'title'), 'string', topoplot_titles{b}, 'FontSize', 14);
            end
            b = b + 1;
        end
        clear a b
            
        % update axis counter
        axis_counter = axis_counter + 1;
    end                                       

    % finalize the summary figure
    figure(fig)  
    subplot(length(CAPSTEP_TEP_default.peaks), 7, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if participant(p) < 10
        subj = ['0' num2str(participant(p))];
    else
        subj = num2str(participant(p));
    end
    figure_name = ['CAPSTEP_amplitude_' subj];
    savefig([folder_figures '\TEP amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\TEP amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
end   
    
% append progressively the output variables to the general MATLAB file
save(output_file, 'CAPSTEP_TEP_peaks', 'CAPSTEP_TEP_default', '-append');
clear p s figure_title k c t e data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak a ...
    col_fig1 col_fig pos finish axis_counter topoplot_titles eoi  

%% 5) save data in a R-compatible table 
if ~exist('CAPSTEP_TEP_values')
    CAPSTEP_TEP_values = table;
end

% export
row_counter = height(CAPSTEP_TEP_values) + 1;
for p = 1:length(participant) 
    for c = 1:length(condition)  
        for t = 1:length(time)
            for k = 1:length(CAPSTEP_TEP_default.peaks) 
                %fill in the table
                CAPSTEP_TEP_values.subject(row_counter) = p;
                CAPSTEP_TEP_values.condition(row_counter) = condition(c);
                CAPSTEP_TEP_values.time(row_counter) = time(t);
                CAPSTEP_TEP_values.peak(row_counter) = CAPSTEP_TEP_default.peaks(k);
                CAPSTEP_TEP_values.amplitude_peak(row_counter) = CAPSTEP_TEP_peaks.amplitude_peak(c, t, p, k);
                CAPSTEP_TEP_values.amplitude_mean(row_counter) = CAPSTEP_TEP_peaks.amplitude_mean(c, t, p, k);
                CAPSTEP_TEP_values.latency(row_counter) = CAPSTEP_TEP_peaks.latency(c, t, p, k);

                % update the counter
                row_counter = row_counter + 1;
            end
        end
    end
end
writetable(CAPSTEP_TEP_values, [folder_results '\CAPSTEP_TEP_values.csv'])

clear c t p k row_counter

%% 6) plot baseline TEPs 
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
if exist('CAPSTEP_TEP_data') ~= 1
    load(output_file, 'CAPSTEP_TEP_peaks', 'CAPSTEP_TEP_default')
end

% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% average baseline data across sessions
data_visual = squeeze(mean(CAPSTEP_TEP_data(:, 1, WO_idx, :, :), [1 3]));

% plot butterfly + highlight Cz
fig = plot_TEP(x, data_visual, time_window, figure_counter);

% name and save figure
figure_name = 'CAPSTEP_TEP_baseline';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])

% update figure counter
figure_counter = figure_counter + 1 ;

clear outliers WO_idx data_visual fig figure_name

%% 7) plot peak values - boxplot
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
if exist('CAPSTEP_TEP_peaks') ~= 1
    load(output_file, 'CAPSTEP_TEP_peaks', 'CAPSTEP_TEP_default')
end

% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% loop through peaks
for k = 1:length(CAPSTEP_TEP_default.peaks)  
    % load data
    data_amplitude = []; data_latency = [];
    for c = 1:length(condition)       
        for t = 1:length(time)
            data_amplitude(c, t, :) = CAPSTEP_TEP_peaks.amplitude_norm(c, t, WO_idx, k);
            data_latency(c, t, :) = CAPSTEP_TEP_peaks.latency(c, t, WO_idx, k);
        end
    end

    % plot amplitude
    fig = plot_box(data_amplitude, 'amplitude', condition, colours, figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('CAPSTEP_TEP_amplitude_%s_WO', CAPSTEP_TEP_default.peaks{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
    
    % plot latency
    fig = plot_box(data_latency, 'latency', condition, colours, figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('CAPSTEP_TEP_latency_%s_WO', CAPSTEP_TEP_default.peaks{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear k c t fig figure_name data_amplitude data_latency  WO_idx outliers

%% 8) plot ratings - lineplot
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
if exist('CAPSTEP_ratings') ~= 1
    load(output_file, 'CAPSTEP_ratings')
end

% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% prepare y, calculate SEM or 95% CI
caps_idx = find(cellfun(@(x) all(x == 1), CAPSTEP_ratings(:, 3)));
ratings_caps = cell2mat(CAPSTEP_ratings(caps_idx, [4:2:30]));
ratings(1, :, :) = ratings_caps(WO_idx, :);
ctrl_idx = find(cellfun(@(x) all(x == 0), CAPSTEP_ratings(:, 3)));
ratings_ctrl = cell2mat(CAPSTEP_ratings(ctrl_idx, [4:2:30]));
ratings(2, :, :) = ratings_ctrl(WO_idx, :);
clear ratings_caps ratings_ctrl
for c = 1:length(condition)
    for t = 1:size(ratings, 3)
        y(c, t) = mean(squeeze(ratings(c, :, t)));
    %     CI(c, t) = std(squeeze(ratings(c, :, t))) / sqrt(size(ratings, 2)) * z;
        CI(c, t) = std(squeeze(ratings(c, :, t))) / sqrt(size(ratings, 2));
    end
end

% prepare x
x = [0.75:0.5:7.25];

% launch the figure
fig = figure(figure_counter);
hold on

% determine limits
xl = [0.25, length(time) + 0.75];
xlim(xl); ylim([-10, 100])

% plot
plot_line(x, y, CI, colours, alpha)

% add parameters
figure_title = sprintf('Intensity ratings\nwithout outliers');
set(gca, 'xtick', 1:length(time), 'xticklabel', time)
set(gca, 'Fontsize', 16)
% title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
xlabel('timepoint'); ylabel(sprintf('intensity\n(NRS +- SEM)'));
set(gca, 'ytick', 0:20:100, 'YColor', [0 0 0])

% change figure size
fig.Position = [500 500 600 250];

% name and save figure
figure_name = sprintf('CAPSTEP_ratings_WO');
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

% update the counter
figure_counter = figure_counter + 1;      
clear t x y caps_idx ctrl_indx CI fig xl P F figure_title figure_name ratings WO_idx outliers

%% 9) plot peak values - lineplot
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
load(output_file, 'CAPSTEP_TEP_peaks', 'CAPSTEP_TEP_default')

% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% prepare x
x_TEP = 1:length(time);

% loop through peaks
for k = 1:length(CAPSTEP_TEP_default.peaks)  
    % prepare y - TEPs
    y_TEP_all = CAPSTEP_TEP_peaks.amplitude_norm(:, :, WO_idx, k);  
    for c = 1:length(condition)        
        for t = 1:length(time)
            y_TEP(c, t) = mean(squeeze(y_TEP_all(c, t, :)));
%             CI_TEP(c, t) = std(squeeze(y_TEP_all(c, t, :))) / sqrt(size(y_TEP_all, 3)) * z;
            CI_TEP(c, t) = std(squeeze(y_TEP_all(c, t, :))) / sqrt(size(y_TEP_all, 3));
        end
    end
    
    % launch the figure
    fig = figure(figure_counter); 
    hold on
    
    % determine x axis properties
    xl = [0.25, length(time) + 0.75];
    xlim(xl);
    xlabel('timepoint');  
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    
    % add no-change line
    line(xl, [0, 0], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.2);
    
    % plot TEP peak values    
    plot_line(x_TEP, y_TEP, CI_TEP, colours, alpha, 'legend', {'capsaicin' 'control'})
    
    % determine y axis properties
    yl_old = get(gca, 'ylim');
    yl_new = yl_old;
%     yl_new(1) = yl_old(1) - round((yl_old(2) - yl_old(1))/4, 1);
    ylim(yl_new)
    set(gca, 'ytick', yl_old(1):0.5:yl_old(2), 'YColor', [0 0 0])
    ylabel('\Delta TEP amplitude (\muV \pm SEM)');
    
    % change figure size
    fig.Position = [500 300 600 400];

    % add other parameters
    set(gca, 'Fontsize', 14)
    figure_title = sprintf('TEP peak amplitude - %s\nwithout outliers', CAPSTEP_TEP_default.peaks{k});
%     title(figure_title, 'FontWeight', 'bold', 'FontSize', 12)

    % name and save figure
    figure_name = sprintf('CAPSTEP_TEP_%s_WO', CAPSTEP_TEP_default.peaks{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update the counter
    figure_counter = figure_counter + 1;    
end
clear k c t x_ratings y_ratings CI_ratings x_TEP y_TEP_all y_TEP CI_TEP caps_idx fig xl yl_old yl_new ...
    figure_title figure_name ratings_caps WO_idx outliers

%% 10) plot peak values + ratings - lineplot
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
load(output_file, 'CAPSTEP_ratings', 'CAPSTEP_TEP_peaks')

% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% prepare y - ratings
caps_idx = find(cellfun(@(x) all(x == 1), CAPSTEP_ratings(:, 3)));
ratings_caps = cell2mat(CAPSTEP_ratings(caps_idx, [4:2:30]));
ratings_caps = ratings_caps(WO_idx, :);
for t = 1:size(ratings_caps, 2)
    y_ratings(t) = mean(ratings_caps(:, t));
%     CI_ratings(t) = std(ratings_caps(:, t)) / sqrt(size(ratings_caps, 1)) * z;
    CI_ratings(t) = std(ratings_caps(:, t)) / sqrt(size(ratings_caps, 1));
end

% prepare x
x_TEP = 1:length(time);
x_ratings = [0.75:0.5:7.25];

% loop through peaks
for k = 1:length(CAPSTEP_TEP_default.peaks)  
    % prepare y - TEPs
    y_TEP_all = CAPSTEP_TEP_peaks.amplitude_norm(:, :, WO_idx, k);  
    for c = 1:length(condition)        
        for t = 1:length(time)
            y_TEP(c, t) = mean(squeeze(y_TEP_all(c, t, :)));
%             CI_TEP(c, t) = std(squeeze(y_TEP_all(c, t, :))) / sqrt(size(y_TEP_all, 3)) * z;
            CI_TEP(c, t) = std(squeeze(y_TEP_all(c, t, :))) / sqrt(size(y_TEP_all, 3));
        end
    end
    
    % launch the figure
    fig = figure(figure_counter); 
    hold on
    
    % determine x axis properties
    xl = [0.25, length(time) + 0.75];
    xlim(xl);
    xlabel('timepoint');  
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    
    % plot ratings
    yyaxis right
    plot_line(x_ratings, y_ratings, CI_ratings, [0 0 0], alpha)
    
    % determine y axis properties
    ylim([-10, 400])
    set(gca, 'ytick', 0:20:100, 'YColor', [0 0 0])
    ylabel('perceived intensity (\pm SEM)');
    
    % add no-change line
    yyaxis left
    line(xl, [0, 0], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.2);
    
    % plot TEP peak values    
    plot_line(x_TEP, y_TEP, CI_TEP, colours, alpha, 'legend', {'capsaicin' 'control'})
    
    % determine y axis properties
    yl_old = get(gca, 'ylim');
    yl_new = yl_old;
    yl_new(1) = yl_old(1) - round((yl_old(2) - yl_old(1))/4, 1);
    ylim(yl_new)
    set(gca, 'ytick', yl_old(1):0.5:yl_old(2), 'YColor', [0 0 0])
    ylabel('\Delta amplitude (\muV \pm SEM)');
    
    % add other parameters
    set(gca, 'Fontsize', 12)
%     figure_title = sprintf('TEP peak amplitude - %s\nwithout outliers', CAPSTEP_TEP_default.peaks{k});
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 12)

    % name and save figure
    figure_name = sprintf('CAPSTEP_TEP+ratings_%s_WO', CAPSTEP_TEP_default.peaks{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update the counter
    figure_counter = figure_counter + 1;    
end
clear k c t x_ratings y_ratings CI_ratings x_TEP y_TEP_all y_TEP CI_TEP caps_idx fig xl yl_old yl_new ...
    figure_title figure_name ratings_caps WO_idx outliers

%% functions
function peak_x = gfp_plot(x, y, time_window, xstep, labeled, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% launch the figure  
plot(x, y)
yl = get(gca, 'ylim');
cla

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data, mark TMS stimulus
plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 5, 'MinPeakProminence', 0.01);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.015
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.220
        idx(a) = false;        
    else
        idx(a) = true;
    end
end
pks = pks(idx); locs = locs(idx);
if length(pks) > max_peaks
    pks = pks(1:max_peaks); 
    locs = locs(1:max_peaks);
end

% calculate peak coordinations
for a = 1:length(locs)
    peak_x(a) = time_window(1) + locs(a)*xstep;
    peak_y(a) = pks(a);
end
peak_y = double(peak_y);

% plot peaks
for a = 1:length(locs)
    plot(peak_x(a), peak_y(a), ...
        'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    line([peak_x(a), peak_x(a)], [yl(1), peak_y(a)], 'LineStyle', ':', 'Color', [1, 0, 0], 'LineWidth', 1.5)
    
    % label with latency (ms)
    if strcmp(labeled, 'on') 
        text(peak_x(a), peak_y(a) + 0.15, sprintf('%1.0fms', peak_x(a)*1000), 'Color', [1 0 0], 'FontSize', 14)
    end
end

% add parameters
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
ylabel('potential (\muV)')
end
function topo_plot(header, data, x_pos, x_start, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - x_start)/header.xstep);
vector = data(1, :, 1, 1, 1, x_visual);

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i=1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i=i+1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
function [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, xstart, polarity)
% identify TOI data
start = ceil(((center-span/2) - xstart)/xstep) + 1;
stop = ceil(((center+span/2) - xstart)/xstep);
x_TOI = start : stop;
y_TOI = data_visual(:, x_TOI);

% calculate number of points to average
x_include_n = ceil((percent/100) * size(y_TOI, 2));

% identify peak values for all datasets  
if polarity > 0
    [y_max, x_max] = max(y_TOI, [], 2);
else
    [y_max, x_max] = min(y_TOI, [], 2);
end
x_peak = start + x_max - 1; 

% identify timepoints that will be included in the mean amplitude
x_TOI_include = [];
for a = 1:size(data_visual, 1)
    x_TOI_include_i = x_max(a);
    t_left = ceil((x_include_n - 1)/2); tx_left = x_max(a) - t_left : x_max(a) - 1;
    t_right = x_include_n - 1 - t_left; tx_right = x_max(a) + 1 : x_max(a) + t_right;
    
    % first look left, check for limit
    n = length(find(tx_left <= 0));
    if n > 0
        tx_left = tx_left(tx_left > 0);
        for b = 1:n
            tx_right(end + 1) = tx_right(end) + 1;
        end
    end
    
    % look right, check for limit
    n = length(find(tx_right > size(y_TOI, 2)));
    if n > 0
        tx_right = tx_right(tx_right <= size(y_TOI, 2));
        for b = 1:n
            tx_left(end + 1) = tx_left(1) - b;
        end
    end
    
    % append approved datapoints
    x_TOI_include_i = sort([x_TOI_include_i tx_left tx_right]);
    x_TOI_include = [x_TOI_include; x_TOI_include_i];    
end
x_include = x_TOI_include + start - 1;

% calculate the mean value
for a = 1:size(data_visual, 1)
    y_mean(a, 1) = mean(data_visual(a, x_include(a, :)));
end

% calculate peak latency
lat_peak = x_peak * xstep + xstart;
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
function fig = plot_TEP(x, data_visual, time_window, figure_counter)   
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of y
    plot(x, data_visual, 'LineWidth', 0.5)
    yl = get(gca, 'ylim'); 
    yl(1) = yl(1) - ((yl(2) - yl(1))*0.05); yl(2) = yl(2) + ((yl(2) - yl(1))*0.05);
    clf, hold on

    % shade interpolated interval 
    rectangle('Position', [-0.005, yl(1) + 0.02, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % loop through channels to plot
    for a = 1:size(data_visual, 1)     
        P(a) = plot(x, data_visual(a, :), 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    end

    % Cz electrode
    P(end+1) =  plot(x, data_visual(18, :), 'Color', [0 0 0], 'LineWidth', 3);

    % TMS stimulus
    line([0, 0], yl, 'Color', [0.8000    0.0549    0.0549], 'LineWidth', 3)

    % other parameters
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)    
    ylim(yl)
    
    % add legend
    lgd = legend(P(end), 'Cz electrode', 'Box', 'off');
    lgd.FontSize = 14;
    lgd.Position = [0.225 0.05 0.4 0.3];
end
function fig = plot_box(data, datatype, condition, col, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % determine x limits
    xl = [0.25, size(data, 2)-0.25];

    % add zero line
    line(xl, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.9)

    % plot data per condition
    for c = 1:length(condition)  
        % subset data
        data_visual = squeeze(data(c, 2:end, :))';

        % determine x positions
        if c == 1
            data_x(c, :) = (1:size(data_visual, 2))-0.17;
        else
            data_x(c, :) = (1:size(data_visual, 2))+0.17;
        end

        % boxplot
        for t = 1:size(data_visual, 2)
            P(c, t) = boxchart(data_x(c, t) * ones(size(data_visual, 1), 1), data_visual(:, t), ...
                'BoxFaceColor', col(c, :), 'BoxWidth', 0.3, 'WhiskerLineColor', col(c, :), 'MarkerColor', col(c, :));
        end       
    end

    % add legend
    lgd = legend(P(:, 1), {'capsaicin' 'control'}, 'Location', 'southeast');
    lgd.FontSize = 14;
    legend('boxoff')

    % y label
    if strcmp(datatype, 'amplitude')
        ylabel('change in amplitude (\muV)')
    elseif strcmp(datatype, 'latency')
        ylabel('change in latency (ms)')
    end

    % other parameters
    xlim(xl)
    xlabel('timepoint post application')
    set(gca, 'FontSize', 14) 
    set(gca, 'layer', 'top');
end
function plot_line(x, y, CI, colours, alpha, varargin)
    % check for varargins
    L = find(strcmpi(varargin, 'legend'));
    if ~isempty(L)
        label = varargin{L + 1};
    end

    % loop through datasets
    for a = 1:size(y, 1)
        % plot the CI
        F(a) = fill([x fliplr(x)],[y(a, :) + CI(a, :) fliplr(y(a, :) - CI(a, :))], ...
            colours(a, :), 'FaceAlpha', alpha, 'linestyle', 'none');

        % plot the data
        P(a) = plot(x, y(a, :))
        P(a).Color = colours(a, :);
        P(a).LineStyle = '-';
        P(a).LineWidth = 1.2;
        P(a).Marker = 'o';
        P(a).MarkerFaceColor = colours(a, :);
        P(a).MarkerSize = 6;
    end

    % add legend if required
    if ~isempty(L)
        lgd = legend(P, label, 'Location', 'northwest');
        lgd.FontSize = 14;
        legend('boxoff')
    end
end










