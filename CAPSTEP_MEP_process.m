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

%% parameters
clear all
clc

% dataset
participant = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '12' '13' '14' '16' '17'};
session = {'ctrl' 'caps'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix_old = 'dc ep EMG'; 
prefix_1 = 'visual';
prefix_2 = 'nozero';
output_name = 'CAPSTEP_MEP';

% filter
baseline = -0.2;
threshold = 15;
allow_sd = 3;

% amplitude
window = [0.015 0.05];

% visualization
figure_counter = 1;
z = 1.96;
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
            for s = 1:header.datasize(1)
                header.events(s).epoch = s;
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
            savefig(fig, [figure_name '.fig'])
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

%% PEAK-TO-PEAK AMPLITUDE
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
                
                % loop through epochs
                for e = 1:size(data, 1)
                    % choose the window
                    x_start = ceil((window(1) - header.xstart)/header.xstep);
                    x_end = ceil((window(2) - header.xstart)/header.xstep);
                    
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

% save output table
save([output_name '.mat'], 'MEP')
clear amplitudes

%% 6) VISUALIZATION per session - not normalized 
cmap = {'summer' 'autumn'};
for s = 1:length(session)
    % choose the data
    data_visual = [];
    for t = 1:length(time)
        rows = (categorical(MEP.session) == session{s} & ...
            categorical(MEP.timepoint) == time{t});
        data_i = MEP.amplitude_zero(rows);
        data_visual = cat(2, data_visual, data_i);
    end
    
    % plot group boxplot
    fig = figure(figure_counter);        
    boxplot(data_visual, 'color', 'k')
    hold on
    
    % prepare colours
    statement = ['colormap ' cmap{s}]; eval(statement)
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
    figure_name = ['CAPSTEP_MEP_' session{s}];
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    set(gca, 'Fontsize', 14)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('timepoint'); ylabel('MEP (\muV)');
    
    % save the figure       
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])
    
    % update the counter
    figure_counter = figure_counter + 1;    
end
clear s t p b fig pline scat statement C figure_title figure_name cmap rows data_i data_visual 

%% 7) VISUALIZATION per session - normalized 
% calculate normalized amplitudes
session = session([2 1]);
amp_norm = [];
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time) 
            if t == 1
                % baseline --> 100 %
                amp_norm(end + 1) = 100; 
            else
                % current MEP
                rows = (categorical(MEP.subject) == participant{p} & ...
                    categorical(MEP.session) == session{s} & ...
                    categorical(MEP.timepoint) == time{t});
                data_i = MEP.amplitude_zero(rows);
                
                % baseline MEP
                rows = (categorical(MEP.subject) == participant{p} & ...
                    categorical(MEP.session) == session{s} & ...
                    categorical(MEP.timepoint) == time{1});
                data_bl = MEP.amplitude_zero(rows);
                
                % normalize in % of baseline
                amp_norm(end + 1) = data_i / data_bl * 100;
            end            
        end
    end
end
session = session([2 1]);
clear p s t rows data_i data_bl

% append the amplitudes to the result table
MEP.amplitude_norm = amp_norm';
save([output_name '.mat'], 'MEP')
clear amp_norm

% plot 
cmap = {'summer' 'autumn'};
for s = 1:length(session)
    % choose the data
    data_visual = [];
    for t = 1:length(time)
        rows = (categorical(MEP.session) == session{s} & ...
            categorical(MEP.timepoint) == time{t});
        data_i = MEP.amplitude_norm(rows);
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
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])
    
    % update the counter
    figure_counter = figure_counter + 1;    
end
clear s t p b fig pline scat statement C figure_title figure_name cmap rows data_i data_visual 

%% 8) VISUALIZATION both sessions together - normalized 
% choose the data
data_visual = [];
for s = 1:length(session)
    for t = 1:length(time)
        rows = (categorical(MEP.session) == session{s} & ...
            categorical(MEP.timepoint) == time{t});
        data_i = MEP.amplitude_norm(rows);
        data_visual = cat(2, data_visual, data_i);
    end
end
data_visual = data_visual([1:8, 10:11, 13:15], :);

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
        CI(t) = std(data_visual(:, (s - 1) * length(time) + t)) / sqrt(length(participant)) * z;
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
figure_title = ['MEP normalized - both sessions'];
figure_name = ['CAPSTEP_MEP_both_WO'];
set(gca, 'xtick', 1:length(time), 'xticklabel', time)
set(gca, 'Fontsize', 14)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 16)
xlabel('timepoint'); ylabel('MEP (% baseline)');

% save the figure       
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;    

clear s t rows data_i data_visual x y perr figure_title figure_name CI 





 


               