%% EMG PROCESSING
% Written by Dominika for the CAPS-TEP project
% 
% ----- Discards epochs with baseline activity ----- 
% 1) Loops through the datset, in one loop:
%       - calculates average RMS of the baseline  
%       - discards all epochs that have baseline RMS larger than 
%         average RMS + <allow_sd> * sd
%    In case that there are discarded epochs, recalculates --> next cycle
% 2) Discards all epochs that at any point depass threshold value 
% 3) Saves the dataset, updates outcome table + saves outcome figure
% 
% ----- Extracts MEP peak-to-peak amplitudes ----- 
% 4) Calculates the amplitude for each epoch
% 5) Identifies zero epochs - p2p amplitude < 2 * threshold
%       - splits in two variables with and without zero epochs  
%       - mean p2p amplitude --> in the 'output' table
%       - saves sorted data as LW dataset --> prefix 'nozero'
% 
% ----- Group visualization ----- 
% 6) Plots group values 


%% parameters
clear all
clc

% dataset
participant = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '12' '13' '14' '16' '17'};
session = {'ctrl' 'caps'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};

% output
prefix_old = '';
prefix_1 = 'clean';
prefix_2 = 'nozero';
output_name = 'CAPS-TEP_MEP';

% filter
baseline = -0.2;
threshold = 15;
allow_sd = 3;

% amplitude
window = [0.015 0.05];

% visualization
figure_counter = 1;
for c = 1:length(session)
    color_scheme(c, :) = uisetcolor([session{c} ' color'])
end
save('color_scheme.mat', 'color_scheme')
clear c 

%% BASELINE ACTIVITY 
% prepare an output table
if exist([output_name '.mat']) == 0
    MEP_results = table;
    MEP_results.subject = zeros(0); MEP_results.session = zeros(0); MEP_results.timepoint = zeros(0);
    MEP_results.discarded = zeros(0); MEP_results.percent = zeros(0); MEP_results.cycles = zeros(0); 
    MEP_results.threshold = zeros(0);
else
    load('MEP_filtered.mat')
end
    
% loop through participants, sessions and timepoints
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)                               
            % load header and dataset 
            load([prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'],'-mat');      
            load([prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.mat']); 

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
                title(['Original dataset - average RMS ' num2str(avg_rms) ' µV'])
                set(gca,'Fontsize',16); ylabel('amplitude (µV)');
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
                set(gca,'Fontsize',16); ylabel('amplitude (µV)');
                hold on
                
                subplot(3, 2, [5 6])
                title(['Discarded epochs: ' num2str(length(discarded_pos))])
                set(gca,'Fontsize',16); xlabel('time (s)'); ylabel('amplitude (µV)');
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
                    set(gca,'Fontsize',14); ylabel('amplitude (µV)');
                    hold on  
                    
                    % plot filtered dataset
                    subplot(3, 2, 5)
                    plot(x, data_visual, 'Color', [0, 0, 0], 'linewidth', 1.5)
                    xlim = get(gca,'xlim');                        
                    title(['RMS + SD: ' num2str(length(header.events)) ' epochs kept, ' num2str(cycle - 1) ' cycles performed'])
                    set(gca,'Fontsize',14); ylabel('amplitude (µV)'); xlabel('time (s)');
                    hold on

                    % add threshold 
                    subplot(3, 2, 5)
                    l(1) = line([-0.2, 0], [threshold, threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                    l(2) = line([-0.2, 0], [-threshold, -threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                    text(xlim(1) + 0.005 , - threshold + 4 ,['threshold = ' num2str(threshold) ' µV'], 'Fontsize', 14, 'color', [0.99, 0.3, 0.2])
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
%                         text(xlim(1) + 0.005 , ylim(2) - 2 ,['Final average RMS: ' num2str(avg_rms) ' µV'], 'Fontsize', 16)
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
            set(gca,'Fontsize',14); ylabel('amplitude (µV)');  
            title(['Subject ' participant{p} ', ' session{s} ', ' time{t} ' - FINAL DATASET : ' num2str(length(header.events) - length(discarded_pos)) ' epochs kept, ' ...
                num2str(length(discarded)) ' discarded - final average RMS ' num2str(final_rms) ' µV'], 'Fontsize', 18)            
            hold on
            
            subplot(3, 2, 6)
            title(['THRESHOLD : ' num2str(length(discarded_pos)) '  epochs discarded'])
            set(gca,'Fontsize',14); xlabel('time (s)');
            hold on
            
            pause(3)
                                            
            %% 3) save outcome variables           
            % number final epochs
            for n = 1:header.datasize(1)
                header.events(n).epoch = n;
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

            MEP_results = [MEP_results; M];

            save([output_name '.mat'], 'MEP_results') 
            
            % save and close the figure
            figure_name = [prefix_1 '_' participant{p} '_' session{s} '_' time{t}];
            savefig(fig, [figure_name '.fig'])
            close(fig)

            % update the counter
            figure_counter = figure_counter + 1;
        end
    end
end
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
    for n = 1:length(session)
        for t = 1:length(time)             
            for s = 1:length(stimulus)
                %% 4) extract p2p amplitude
                % load data and header
                load([prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG.mat']);
                load([prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG.lw6'], '-mat');
                
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
                A.subject = participant(p); A.session = session(n); A.time = time(t); A.stimulus = stimulus(s);
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
end
clear A data header data_nozero x_start x_end 
clear p n t s e

% make sure participants are ordered in result tables
MEP_results = sortrows(MEP_results, [1:3]);
amplitudes = sortrows(amplitudes, [1:3]);

% append results
MEP_results.epochs_zero = amplitudes.zero;
MEP_results.amplitude_zero = amplitudes.amp_zero;
MEP_results.epochs_nozero = amplitudes.nozero;
MEP_results.amplitude_nozero = amplitudes.amp_nozero;

% save output table
save([output_name '.mat'], 'MEP_results')
clear amplitudes

               