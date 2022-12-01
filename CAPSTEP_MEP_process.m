%% EMG DATA PROCESSING, MEP EXTRACTION
% 
% ----- Discards MEP epochs with baseline activity ----- 
% 1) Loops through the datset, in one loop:
%       - calculates average RMS of the baseline  
%       - discards all epochs that have baseline RMS larger than 
%         average RMS + <allow_sd> * sd
%    In case that there are discarded epochs, recalculates --> next cycle
% 2) Discards all epochs that at any point depass <threshold> value 
% 3) Saves the dataset, updates MEP result table + saves outcome figure

% ----- Extracts peak-to-peak amplitude ----- 
% 4) Calculates the amplitude for each epoch
% 5) Identifies zero epochs - p2p amplitude < 2 * threshold
%       - splits in two variables with and without zero epochs  
%       - mean p2p amplitude --> in the MEP result table
%       - saves sorted data as LW dataset --> prefix 'nozero'
% 
% ----- Group visualization ----- 
% 6) Plots each session separately - absolute MEP amplitude
%       - boxplot + scatterplot
%       - visualizes each subject individually to see the trend
% 7) Plots each session separately - normalized MEP amplitude
%       - calculates amplitudes normalized as % of baseline, appends the
%       values to MEP result table
%       - line + scatterplot, individual subjects
% 8) Plots both sessions together - normalized MEP amplitude
%       - dotplot
% 9) Plots both sessions together (normalized) + adds ratings
%       - lineplot + SEM shading
% 10) Plots average MEP timecourse
%       - at the baseline
%       - baseline vs. t6

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
output_file = [folder_results '\CAPSTEP_output'];

% dataset
participant = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '12' '13' '14' '16' '17' '18' '20' '22' '23' '24'};
session = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix_old = 'dc ep EMG'; 
prefix_1 = 'visual';
prefix_2 = 'nozero';
output_name = 'CAPSTEP_MEP_2';

% filter
baseline = -0.2;
threshold = 15;
allow_sd = 3;

% amplitude
window = [0.015 0.05];

% visualization
figure_counter = 1;
z = 1.96;
alpha = 0.15;

load([folder_git '\colours.mat'])
if ~exist('colormap.mat')
    for c = 1:length(session)
        col(c, :) = uisetcolor(['Choose colour for ' session{c} ' data:']);
    end
    save('colormap.mat', 'col')
    clear c
else
    load('colormap.mat')
end

%% BASELINE ACTIVITY
% p = 1; s = 2; t = 6;    
% prepare an output table
if exist([output_name '.mat']) == 0
    MEP = table;
    MEP.subject = zeros(0); MEP.session = zeros(0); MEP.timepoint = zeros(0);
    MEP.discarded = zeros(0); MEP.percent = zeros(0); MEP.cycles = zeros(0); 
    MEP.threshold = zeros(0);
else
    load([output_name '.mat'])
end
    
% loop through participants, sessions and timepoints
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)                               
            % load header and dataset 
            load([prefix_old ' ' char(participant(p)) ' ' char(session(s)) ' ' char(time(t)) '.lw6'],'-mat');      % load the header
            load([prefix_old ' ' char(participant(p)) ' ' char(session(s)) ' ' char(time(t)) '.mat']); 

            % prepare baseline data for visualization, detrend
            data_visual = squeeze(data(:, 1, 1, 1, 1, 1 : ceil(-baseline / header.xstep)));
            for a = 1:size(data_visual, 1)
                data_visual(a, :) = detrend(data_visual(a, :));
            end
            data_visual_orig = data_visual;
            
            % prepare x 
            x = baseline : header.xstep : 0; 
            
            % outcome vector
            discarded = [];

            %% 1) run automatic RMS + SD removal 
            go = 1; cycle = 1; 
            while go 
                discarded_pos = [];

                % calculate average RMS values
                avg_rms = mean(rms(data_visual'));
                avg_sd = std(rms(data_visual'));
                cutoff = avg_rms + allow_sd * avg_sd;

                % plot original dataset
                fig = figure(figure_counter); 
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                subplot(3, 2, [1 2])
                plot(x, data_visual, 'color', [0.45, 0.45, 0.45], 'linewidth', 1.5)                
                sgtitle(['Subject ' participant{p} ' - ' session{s} ' - ' time{t} ' : CYCLE ' num2str(cycle)], 'fontsize', 18)
                title(['Original dataset - average RMS ' num2str(avg_rms) ' �V'])
                set(gca,'Fontsize',16); ylabel('amplitude (�V)');
                hold on
                              
                % loop through trials
                for e = 1:length(header.events)
                    % check if the rms of current event fits into the limis 
                        current_rms = rms(data_visual(e, :)');
                        if current_rms > cutoff
                            discarded_pos = [discarded_pos e];
                            discarded = [discarded header.events(e).epoch];
                            ditch = true;
                        else
                            ditch = false;
                        end                        
                       
                    % plot the event to the appropriate axes                
                    if ditch   
                        subplot(3, 2, [5 6])
                        plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
                        hold on
                    else
                        subplot(3, 2, [3 4])
                        plot(x, data_visual(e, :), 'Color', [0, 0, 0], 'LineWidth', 1.5)
                        hold on
                    end
                end

                % add parameters to axes
                subplot(3, 2, [3 4])
                title(['Kept epochs: ' num2str(length(header.events) - length(discarded_pos))])
                set(gca,'Fontsize',16); ylabel('amplitude (�V)');
                hold on
                
                subplot(3, 2, [5 6])
                title(['Discarded epochs: ' num2str(length(discarded_pos))])
                set(gca,'Fontsize',16); xlabel('time (s)'); ylabel('amplitude (�V)');
                hold on

                % ------------------ cycle automatically to 0 discarded epochs ------------------
                if ~isempty(discarded_pos)
                    % remove indicated epochs from data, update header
                    data(discarded_pos, :, :, :, :, :) = [];
                    header.datasize(1) = header.datasize(1) - length(discarded_pos);
                    header.events(discarded_pos)= [];

                    % remove indicated epochs from visual dataset for
                    % future filtration cycles
                    data_visual(discarded_pos , :) = [];
                    
                    % continue the cycle
                    pause(1); clf;
                    cycle = cycle + 1;  
                else
                    % close current figure
                    pause(1); close(fig);   

                    % plot original dataset
                    fig = figure(figure_counter); 
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                    subplot(3, 2, 3)
                    plot(x, data_visual_orig, 'color', [0, 0, 0], 'linewidth', 1.5)                       
                    title(['Original dataset: ' num2str(size(data_visual_orig, 1)) ' epochs'])
                    set(gca,'Fontsize',14); ylabel('amplitude (�V)');
                    hold on  
                    
                    % plot filtered dataset
                    subplot(3, 2, 5)
                    plot(x, data_visual, 'Color', [0, 0, 0], 'linewidth', 1.5)
                    xlim = get(gca,'xlim');                        
                    title(['RMS + SD: ' num2str(length(header.events)) ' epochs kept, ' num2str(cycle - 1) ' cycles performed'])
                    set(gca,'Fontsize',14); ylabel('amplitude (�V)'); xlabel('time (s)');
                    hold on

                    % add threshold 
                    subplot(3, 2, 5)
                    l(1) = line([-0.2, 0], [threshold, threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                    l(2) = line([-0.2, 0], [-threshold, -threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                    text(xlim(1) + 0.005 , - threshold + 4 ,['threshold = ' num2str(threshold) ' �V'], 'Fontsize', 14, 'color', [0.99, 0.3, 0.2])
                    hold on
                    
                    % plot discarded epochs - RMS + SD
                    if ~isempty(discarded)
                        subplot(3, 2, 4)
                        plot(x, data_visual_orig(discarded, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5) 
                        title(['RMS + SD: ' num2str(length(discarded)) ' epochs discarded'])
                        set(gca,'Fontsize',14)
                        hold on
                    else
                        subplot(3, 2, 4)
                        title('No epoch discarded'); set(gca,'Fontsize',14)
                        hold on
                    end                     

                    % exit the while loop
                    go = 0;
                end
                
                % --------------------- choose final cycles manually ---------------------
%                     answer = questdlg('Do you want to repeat the cycle?', 'Cycles',...
%                     'Yes, keep them comin', 'No, save and continue to next dataset', 'No, save and continue to next dataset');
%                     switch answer
%                         case 'Yes, keep them comin'
%                             cycle = cycle + 1; 
%                             clf
%                         case 'No, save and continue to next dataset' 
%                             % close currentfigure
%                             close(fig)     
%                             
%                          % plot the final figure and wait 5s
%                         fig = figure(figure_counter); 
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%                         subplot(3, 2, [1 2])
%                         plot(x, data_visual_orig)
%                         h = suptitle(['Subject ' participant{p} ' - ' session{s} ' - ' time{t} ' : ' num2str(cycle - 1) ' cycles performed'])
%                         set(h,'FontSize',20)
%                         title(['Original dataset - ' num2str(size(data_visual_orig, 1)) ' epochs'])
%                         set(gca,'Fontsize',16)
%                         hold on  
% 
%                         subplot(3, 2, [3 4])
%                         plot(x, data_visual, 'Color', [0, 0, 0], 'LineWidth', 1)
%                         ylim = get(gca,'ylim'); xlim = get(gca,'xlim');                        
%                         text(xlim(1) + 0.005 , ylim(2) - 2 ,['Final average RMS: ' num2str(avg_rms) ' �V'], 'Fontsize', 16)
%                         title(['Kept epochs: ' num2str(length(header.events))])
%                         set(gca,'Fontsize',16)
%                         hold on
% 
%                         if ~isempty(discarded)
%                             subplot(3, 2, [5 6])
%                             plot(x, data_visual_orig(discarded, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
%                             title(['Discarded epochs: ' num2str(length(discarded))])    
%                             set(gca,'Fontsize',16)
%                             hold off
%                         end 
%                             
%                             pause(5)
%                             
%                             % save and close the figure
%                             figure_name = ['filtered_' method{mode} '_' participant{p} '_' session{s} '_' time{t}]
%                             savefig(fig, [figure_name '.fig'])
%                             close(fig)
%                             
%                             % update the counter
%                             figure_counter = figure_counter+ 1;
%                             
%                             % exit the while loop
%                             go = 0;
%                     end
%                     % -----------------------------------------------------------------------
            end
            
            %% 2) remove epochs that depass the threshold
            % loop through left trials
            for e = 1:length(header.events)
                % check if tthe maximum value across baseline datapoints
                % fits under the threshold
                    current_max = max(abs(data_visual(e, :)));
                    if current_max > threshold
                        discarded_pos = [discarded_pos e];
                        discarded = [discarded header.events(e).epoch];
                        ditch = true;
                    else
                        ditch = false;
                    end                        

                % plot the event to the appropriate axes                
                if ~ditch   
                    subplot(3, 2, [1 2])
                    plot(x, data_visual(e, :), 'Color', [0, 0.45, 0.74], 'LineWidth', 1.5)
                    hold on
                else
                    subplot(3, 2, 6)
                    plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
                    hold on
                end
            end
                      
            % remove indicated epochs from visual dataset
            data_visual(discarded_pos , :) = [];

            % add parameters to axes
            final_rms = mean(rms(data_visual'));
            subplot(3, 2, [1 2])
            set(gca,'Fontsize',14); ylabel('amplitude (�V)');  
            title(['Subject ' participant{p} ', ' session{s} ', ' time{t} ' - FINAL DATASET : ' num2str(length(header.events) - length(discarded_pos)) ' epochs kept, ' ...
                num2str(length(discarded)) ' discarded - final average RMS ' num2str(final_rms) ' �V'], 'Fontsize', 18)            
            hold on
            
            subplot(3, 2, 6)
            title(['THRESHOLD : ' num2str(length(discarded_pos)) '  epochs discarded'])
            set(gca,'Fontsize',14); xlabel('time (s)');
            hold on
            
            pause(3)
                                            
            %% 3) save outcome variables           
            % number final epochs
            for h = 1:header.datasize(1)
                header.events(h).epoch = h;
            end
            
            % modify and save header
            header.datasize(1) = header.datasize(1) - length(discarded_pos);
            header.events(discarded_pos)= [];
            header.events = header.events';
            header.name = [prefix_1 ' ' prefix_old ' ' participant{p} ' ' session{s} ' ' time{t}];
            save([header.name '.lw6'], 'header')         
           
            % modify and save data
            data(discarded_pos, :, :, :, :, :) = [];
            save([header.name '.mat'], 'data')

            % fill in the outcome table and save
            M = table;
            M.subject = participant(p); M.session = session(s); M.timepoint = time(t); 
            M.discarded = {sort(discarded)}; 
            M.percent = round(100 - (size(data_visual, 1) / size(data_visual_orig, 1)) * 100, 2); 
            M.cycles = cycle - 1; 
            M.threshold = threshold;

            MEP = [MEP; M];

            save([output_name '.mat'], 'MEP') 
            
            % save and close the figure
            figure_name = [prefix_1 '_' participant{p} '_' session{s} '_' time{t}];
            savefig([folder_figures '\' figure_name '.fig'])
            saveas(fig, [folder_figures '\' figure_name '.png'])
            close(fig)

            % update the counter
            figure_counter = figure_counter + 1;
        end
    end
end

% make sure participants are ordered 
MEP = sortrows(MEP, [1:3]);
save([output_name '.mat'], 'MEP') 

xstep = header.xstep;
clear avg_rms avg_sd current_max current_rms cutoff cycle data data_visual data_visual_orig dataset_name ...
    discarded discarded_pos ditch fig figure_name final_rms go header M xlim l x
clear a e i n p s t 

%% PEAK-TO-PEAK AMPLITUDE - subject average
% prepare a table for amplitude results
amplitudes = table; 
amplitudes.subject = zeros(0); amplitudes.session = zeros(0); amplitudes.timepoint = zeros(0);
amplitudes.zero = zeros(0); amplitudes.amp_zero = zeros(0); amplitudes.nozero = zeros(0); amplitudes.amp_nozero = zeros(0);

% loop through datasets 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)             
                %% 4) extract p2p amplitude
                % load data and header
                load([prefix_1 ' ' prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.mat']);
                load([prefix_1 ' ' prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'], '-mat');
                
                % choose the window
                x_start = ceil((window(1) - header.xstart)/header.xstep);
                x_end = ceil((window(2) - header.xstart)/header.xstep);
                
                % loop through epochs
                for e = 1:size(data, 1)                    
                    % identify extremes
                    y_max(e) = max(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                    y_min(e) = min(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                    
                    % calculate amplitude 
                    amp(e) = y_max(e) - y_min(e); 
                end                
                
                %% 5) identify zero response epochs                 
                % loop through epochs
                for e = 1:length(amp)
                    if amp(e) > 2 * threshold
                        zero_index(e) = true; 
                    else
                        zero_index(e) = false; 
                    end
                end
                
                % create new dataset with no zero epochs
                data_nozero = data(zero_index, :, :, :, :, :); 
                
                % extract and append mean values
                A = table;
                A.subject = participant(p); A.session = session(s); A.timepoint = time(t); 
                A.zero = size(data, 1); A.amp_zero = mean(amp); 
                A.nozero = size(data_nozero, 1); A.amp_nozero = mean(amp(zero_index)); 
                
                amplitudes = [amplitudes; A];
                
                % adjust and save header
                header.name = [prefix_2 ' ' header.name];
                header.datasize(1) = size(data_nozero, 1);
                header.events = header.events(1:header.datasize(1));
                save([header.name '.lw6'], 'header')
                
                % save new dataset 
                data = data_nozero;
                save([header.name '.mat'], 'data')
                
                clear y_max y_min amp zero_index            
        end 
    end
end

clear A data header data_nozero x_start x_end 
clear e p t s e

% make sure participants are ordered 
amplitudes = sortrows(amplitudes, [1:3]);

% append results
MEP.epochs_zero = amplitudes.zero;
MEP.amplitude_zero = amplitudes.amp_zero;
MEP.epochs_nozero = amplitudes.nozero;
MEP.amplitude_nozero = amplitudes.amp_nozero;

% append output table to the general MATLAB file
CAPSTEP_MEP = MEP;
save(output_file, 'CAPSTEP_MEP', '-append');
clear amplitudes

%% PEAK FEATURES - single trial
% launch the table 
CAPSTEP_MEP_singletrial = table; 

% loop through datasets 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time) 
            % load data and header
            load([folder_input '\' prefix_1 ' ' prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.mat']);
            load([folder_input '\' prefix_1 ' ' prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'], '-mat');
            
            % choose the window
            x_start = ceil((window(1) - header.xstart)/header.xstep);
            x_end = ceil((window(2) - header.xstart)/header.xstep);

            % extract peak features
            for e = 1:size(data, 1)
                % identify extremes
                y_max = max(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                y_min = min(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                
                % determine peak latency
                x_max = find(data(e, 1, 1, 1, 1, :) == y_max) * header.xstep + header.xstart;
                x_min = find(data(e, 1, 1, 1, 1, :) == y_min) * header.xstep + header.xstart;
                if x_max < x_min
                    lat = round(x_max * 1000, 2);
                else
                    lat = round(x_min * 1000, 2);
                end
                
                % calculate peak amplitude 
                amp = y_max - y_min; 
                
                % append the line to the output table 
                if size(CAPSTEP_MEP_singletrial, 1) == 0
                    CAPSTEP_MEP_singletrial.subject = participant(p); 
                    CAPSTEP_MEP_singletrial.session = session(s); 
                    CAPSTEP_MEP_singletrial.timepoint = time(t); 
                    CAPSTEP_MEP_singletrial.trial = e; 
                    CAPSTEP_MEP_singletrial.latency = lat;
                    CAPSTEP_MEP_singletrial.amplitude = amp;
                else
                    CAPSTEP_MEP_singletrial.subject(end + 1) = participant(p); 
                    CAPSTEP_MEP_singletrial.session(end) = session(s); 
                    CAPSTEP_MEP_singletrial.timepoint(end) = time(t); 
                    CAPSTEP_MEP_singletrial.trial(end) = e; 
                    CAPSTEP_MEP_singletrial.latency(end) = lat;
                    CAPSTEP_MEP_singletrial.amplitude(end) = amp;
                end
            end             
        end
    end
end

% append output table to the general MATLAB file
save(output_file, 'CAPSTEP_MEP_singletrial', '-append');

clear p s t e data header x_start x_end y_max y_min x_max x_min lat amp

%% 6) VISUALIZATION per session - not normalized 
% load data
load(output_file, 'CAPSTEP_MEP')

% plot
cmap = {'autumn' 'winter'};
for s = 1:length(session)
    % choose the data
    data_visual = [];
    for t = 1:length(time)
        rows = (categorical(CAPSTEP_MEP.session) == session{s} & ...
            categorical(CAPSTEP_MEP.timepoint) == time{t});
        data_i = CAPSTEP_MEP.amplitude_zero(rows);
        data_visual = cat(2, data_visual, data_i);
    end
    
    % plot group boxplot
    fig = figure(figure_counter);        
    boxplot(data_visual, 'color', 'k')
    hold on
    
    % prepare colours
    statement = ['colormap ' cmap{s}]; 
    eval(statement)
    C = colormap; 
    C = C(1 : floor(length(C)/length(participant)) : length(C), :);
    
    % plot the lines 
    for p = 1:length(participant)
        pline(p) = plot(1:length(time), data_visual(p, :), '-o',...
            'Color', C(p, :), ...
            'MarkerSize', 10,...
            'MarkerEdge', 'none');
        hold on
    end

    % plot the markers
    for b = 1:length(time)
        scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
            25, C(1:length(participant), :), 'filled');
        hold on
    end
    
    % add parameters
    figure_title = ['MEP amplitude - ' session{s}];
    figure_name = ['CAPSTEP_MEP_' session{s} '_raw'];
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    set(gca, 'Fontsize', 14)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('timepoint'); ylabel('MEP (\muV)');
    
    % save the figure       
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])
    
    % update the counter
    figure_counter = figure_counter + 1;    
end
clear s t p b fig pline scat statement C figure_title figure_name cmap rows data_i data_visual 

%% 7) VISUALIZATION per session - normalized 
% calculate normalized amplitudes
amp_norm = [];
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time) 
            if t == 1
                % baseline --> 100 %
                amp_norm(end + 1) = 100; 
            else
                % current MEP
                rows = (categorical(CAPSTEP_MEP.subject) == participant{p} & ...
                    categorical(CAPSTEP_MEP.session) == session{s} & ...
                    categorical(CAPSTEP_MEP.timepoint) == time{t});
                data_i = CAPSTEP_MEP.amplitude_zero(rows);
                
                % baseline MEP
                rows = (categorical(CAPSTEP_MEP.subject) == participant{p} & ...
                    categorical(CAPSTEP_MEP.session) == session{s} & ...
                    categorical(CAPSTEP_MEP.timepoint) == time{1});
                data_bl = CAPSTEP_MEP.amplitude_zero(rows);
                
                % normalize in % of baseline
                amp_norm(end + 1) = data_i / data_bl * 100;
            end            
        end
    end
end
clear p s t rows data_i data_bl

% append the amplitudes to the result table
CAPSTEP_MEP.amplitude_norm = amp_norm';
save(output_file, 'CAPSTEP_MEP', '-append');
clear amp_norm

% plot 
cmap = {'autumn' 'winter'};
for s = 1:length(session)
    % choose the data
    data_visual = [];
    for t = 1:length(time)
        rows = (categorical(CAPSTEP_MEP.session) == session{s} & ...
            categorical(CAPSTEP_MEP.timepoint) == time{t});
        data_i = CAPSTEP_MEP.amplitude_norm(rows);
        data_visual = cat(2, data_visual, data_i);
    end       

    % prepare colours
    statement = ['colormap ' cmap{s}]; eval(statement)
    C = colormap; 
    C = C(1 : floor(length(C)/length(participant)) : length(C), :);
    
    % plot the lines 
    fig = figure(figure_counter);   
    hold on
    for p = 1:length(participant)
        pline(p) = plot(1:length(time), data_visual(p, :), '-o',...
            'Color', C(p, :), ...
            'MarkerSize', 10,...
            'MarkerEdge', 'none');
    end

    % plot the markers
    for b = 1:length(time)
        scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
            50, C(1:length(participant), :), 'filled');
    end
    
    % add no-change line
    xlim([0.75 7.25])
    line([0.75 7.25], [100, 100], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);
    
    % add parameters
    figure_title = ['MEP normalized - ' session{s}];
    figure_name = ['CAPSTEP_MEP_' session{s} '_norm'];
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    set(gca, 'Fontsize', 14)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('timepoint'); ylabel('MEP (% baseline)');
    
    % save the figure       
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])
    
    % update the counter
    figure_counter = figure_counter + 1;    
end
clear s t p b fig pline scat statement C figure_title figure_name cmap rows data_i data_visual 

%% 8A) VISUALIZATION both sessions together - normalized 
% choose the data
data_visual = [];
for s = 1:length(session)
    for t = 1:length(time)
        rows = (categorical(CAPSTEP_MEP.session) == session{s} & ...
            categorical(CAPSTEP_MEP.timepoint) == time{t});
        data_i = CAPSTEP_MEP.amplitude_norm(rows);
        data_visual = cat(2, data_visual, data_i);
    end
end
data_visual = data_visual([1:8, 10:11, 13:18, 20], :);

% launch the figure
fig = figure(figure_counter); 
hold on

% add no-change line
xlim([0.75 7.25])
line([0.75 7.25], [100, 100], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 2);

% plot data with errorbars
x = 1:length(time);
for s = 1:length(session)
    % calculate the data and 95% CI
    for t = 1:length(time)
        y(t) = mean(data_visual(:, (s - 1) * length(time) + t));
        CI(t) = std(data_visual(:, (s - 1) * length(time) + t)) / sqrt(length(participant([1:8, 10:11, 13:18, 20]))) * z;
    end
    
    % plot
    perr(s) = errorbar(x, y, CI)
    
    % add parameters
    perr(s).Color = col(s, :);
    perr(s).LineWidth = 1.2;
    perr(s).Marker = 'o';
    perr(s).MarkerFaceColor = col(s, :);
    perr(s).MarkerSize = 10;
end

% add parameters
figure_title = sprintf('MEP normalized - both sessions\nwithout outliers');
figure_name = ['CAPSTEP_MEP_both_WO'];
set(gca, 'xtick', 1:length(time), 'xticklabel', time)
set(gca, 'Fontsize', 14)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
xlabel('timepoint'); ylabel('MEP (% baseline)');

% save the figure       
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;    

clear s t rows data_i data_visual x y perr figure_title figure_name CI 

%% 8B) VISUALIZATION both sessions together - normalized 
% ----- section input -----
outliers = [9, 13, 23];
% -------------------------
% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% prepare x
x_MEP = 1:length(time);

% prepare y 
load(output_file, 'CAPSTEP_MEP')
for s = 1:length(session)        
    for t = 1:length(time)
        rows = (categorical(CAPSTEP_MEP.session) == session{s} & ...
            categorical(CAPSTEP_MEP.timepoint) == time{t});
        data_i = CAPSTEP_MEP.amplitude_norm(rows);       
        data_i = data_i(WO_idx);
        y_MEP(s, t) = mean(data_i);
        CI_MEP(s, t) = std(data_i) / sqrt(size(data_i, 1));
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
line(xl, [100, 100], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.2);
    
% plot TEP peak values    
plot_line(x_MEP, y_MEP, CI_MEP, colours, alpha, 'legend', {'capsaicin' 'control'})
    
% determine y axis properties
yl_old = get(gca, 'ylim');
yl_new = yl_old;
%     yl_new(1) = yl_old(1) - round((yl_old(2) - yl_old(1))/4, 1);
ylim(yl_new)
set(gca, 'ytick', yl_old(1):10:yl_old(2), 'YColor', [0 0 0])
ylabel(sprintf('relative MEP amplitude\n(%% baseline +- SEM)'));

% change figure size
fig.Position = [500 300 600 400];

% add other parameters
set(gca, 'Fontsize', 14)
figure_title = 'MEP amplitude - %s\nwithout outliers';
% title(figure_title, 'FontWeight', 'bold', 'FontSize', 14)

% name and save figure
figure_name = 'CAPSTEP_MEP_WO';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

% update the counter
figure_counter = figure_counter + 1;    

clear k c t x_ratings y_ratings CI_ratings x_TEP y_TEP_all y_TEP CI_TEP caps_idx fig xl yl_old yl_new ...
    figure_title figure_name ratings_caps WO_idx outliers


%% 9) VISUALIZATION both sessions together (normalized) + ratings 
% ----- section input -----
participant = [1:10, 12:14, 16:18, 20, 22:24];
outliers = [9, 13, 23];
% -------------------------
% load input if necessary
load(output_file, 'CAPSTEP_ratings', 'CAPSTEP_MEP')

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

% prepare y - MEPs
for s = 1:length(session)        
    for t = 1:length(time)
        rows = (categorical(CAPSTEP_MEP.session) == session{s} & ...
            categorical(CAPSTEP_MEP.timepoint) == time{t});
        data_i = CAPSTEP_MEP.amplitude_norm(rows);       
        data_i = data_i(WO_idx);
        y_MEP(s, t) = mean(data_i);
        CI_MEP(s, t) = std(data_i) / sqrt(size(data_i, 1));
    end
end

% prepare x
x_MEP = 1:length(time);
x_ratings = [0.75:0.5:7.25];
    
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
line(xl, [100, 100], 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1.2);

% plot TEP peak values    
plot_line(x_MEP, y_MEP, CI_MEP, colours, alpha, 'legend', {'capsaicin' 'control'})
    
% determine y axis properties
yl_old = get(gca, 'ylim');
yl_new = yl_old;
yl_new(1) = yl_old(1) - round((yl_old(2) - yl_old(1))/4, 1);
ylim(yl_new)
set(gca, 'ytick', yl_old(1):10:yl_old(2), 'YColor', [0 0 0])
ylabel('\Delta amplitude (% baseline \pm SEM)');
    
% add other parameters
set(gca, 'Fontsize', 12)
figure_title = sprintf('MEP amplitude\nwithout outliers');
title(figure_title, 'FontWeight', 'bold', 'FontSize', 12)

% name and save figure
figure_name = sprintf('CAPSTEP_MEP+ratings_WO');
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

% update the counter
figure_counter = figure_counter + 1;    

clear k s t x_ratings y_ratings CI_ratings x_MEP rows data_i y_MEP CI_MEP caps_idx fig xl yl_old yl_new ...
    figure_title figure_name ratings_caps WO_idx outliers

%% 10) EXTRACT INDIVIDUAL MEAN DATA
% ----- section input -----
participant = [1:10, 12:14, 16:18, 20, 22:24];
time_window = [-0.2 0.1];
% -------------------------
% load random header
load([folder_input '\' prefix_1 ' ' prefix_old ' 01 caps t1.lw6'],'-mat');

% define visualization parameters
x = [time_window(1):header.xstep:time_window(2)];
x_start = floor((time_window(1) - header.xstart)/header.xstep) + 1;
x_end = floor((time_window(2) - header.xstart)/header.xstep) + 1;

% load data, average across epochs
for s = 1:length(session) 
    for t = 1:length(time) 
        for p = 1:length(participant) 
            % define participant
            if participant(p) < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end
            
            % load dataset
            load([folder_input '\' prefix_1 ' ' prefix_old ' ' subj ' ' session{s} ' ' time{t} '.mat'])

            % append averaged data in the matrix
            data_mean = squeeze(mean(data, 1));
            CAPSTEP_MEP_data(s, t, p, :) = data_mean(x_start:x_end);        
        end
    end
end

% save to the global MATLAB file
save(output_file, 'CAPSTEP_MEP_data', '-append');
clear s t p subj header data data_mean x x_start x_end

%% 11) VISUALIZATION MEP timecourse
% ----- section input -----
participant = [1:10, 12:14, 16:18, 20, 22:24];
outliers = [9, 13, 23];
time_window = [-0.025 0.075];
% -------------------------
% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% load random header
load([folder_input '\' prefix_1 ' ' prefix_old ' 01 caps t1.lw6'],'-mat');

% define visualization parameters
x = [time_window(1):header.xstep:time_window(2)];
x_start = floor((time_window(1) - header.xstart)/header.xstep);
x_end = floor((time_window(2) - header.xstart)/header.xstep);

% load data for visualization
for s = 1:length(session) 
    for t = 1:length(time) 
        for p = 1:length(participant) 
            % define participant
            if participant(p) < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end
            
            % load dataset
            load([folder_input '\' prefix_1 ' ' prefix_old ' ' subj ' ' session{s} ' ' time{t} '.mat'])

            % append averaged data in the matrix
            data_mean = squeeze(mean(data, 1));
            data_visual_i(s, t, p, :) = data_mean(x_start:x_end);        
        end
    end
end
data_visual_i = data_visual_i(:, :, WO_idx, :);
clear data

% calculate group average values
data_visual = squeeze(mean(data_visual_i, 3));
CI_visual = squeeze(std(data_visual_i, 0, 3)/sqrt(size(data_visual_i, 3))); % * z;

% plot baseline
data(1, :) = squeeze(mean(data_visual(:, 1, :), 1));
CI(1, :) = squeeze(mean(CI_visual(:, 1, :), 1));
fig = plot_timeseries(x, data, CI, figure_counter, time_window, [0 0 0], alpha)
clear data CI

% name and save figure
figure_name = 'CAPSTEP_MEP_baseline';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])

% update figure counter
figure_counter = figure_counter + 1;

% plot MEP baseline vs. t6 - per session
for s = 1:length(session)
    % get the data
    a = 1;
    for t = [1 4]
        data(a, :) = squeeze(data_visual(s, t, :));
        CI(a, :) = squeeze(CI_visual(s, t, :));
        a = a + 1;
    end
    fig = plot_timeseries(x, data, CI, figure_counter, time_window, cat(1, [0 0 0], colours(s, :)), alpha, 'legend', {'baseline' 't3'})

    % name and save figure
    figure_name = sprintf('CAPSTEP_MEP_change_%s', session{s});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'])

    % update figure counter
    figure_counter = figure_counter + 1;
end

clear s t a outliers WO_idx x x_start x_end data_visual_i data_visual CI_visual data CI fig figure_name

%% FUNCTIONS
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
function fig = plot_timeseries(x, y, CI, figure_counter, time_window, colours, alpha, varargin)
    % check for varargins
    L = find(strcmpi(varargin, 'legend'));
    if ~isempty(L)
        label = varargin{L + 1};
    end
    
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    plot(x, y + CI, 'b:', 'LineWidth', 0.5)
    plot(x, y - CI, 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim'); 
    yl(1) = yl(1) - ((yl(2) - yl(1))*0.05); yl(2) = yl(2) + ((yl(2) - yl(1))*0.05);
    clf, hold on

    % loop through datasets to plot
    for a = 1:size(y, 1)        
        P(a) = plot(x, y(a, :), 'Color', colours(a, :), 'LineWidth', 2.5);
        F(a) = fill([x fliplr(x)],[y(a, :) + CI(a, :) fliplr(y(a, :) - CI(a, :))], ...
            colours(a, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    end

    % mark TMS stimulus
    line([0, 0], yl, 'Color', [0.8000    0.0549    0.0549], 'LineWidth', 3)

    % set other parameters
    xlim(time_window)
    ylim(yl)
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)

    % add legend if required
    if ~isempty(L)
        lgd = legend(P, label, 'Location', 'northeast');
        lgd.FontSize = 12;
        legend('boxoff')
    end
end


               