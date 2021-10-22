%% cut the data
participant = {'12' '13'};
treatment = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix = 'crop bandpass notch prefilt prea dc reref ds art-sup ep dc EEG'; 

% a = 1; b = 1; c = 1;
for a = 1:length(participant)
    for b = 1:length(treatment)
        for c = 1:length(time)
            % deal with the dataset
            load([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.mat']) 
            data = data(:, :, :, :, :, [1:end-1]);
            save([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.mat'], 'data')
            
            % deal with the header
            load([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.lw6'], '-mat') 
            header.datasize(6) = size(data,6);          
            save([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.lw6'], 'header')
        end 
    end
end
clear a b c 

%% rename events
% parameters
participant = {'12' '13'};
treatment = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix = 'crop bandpass notch prefilt prea dc reref ds art-sup ep dc EEG'; 
eventcode = 'Stimulation';

% a = 1; b = 1; c = 1;
for a = 1:length(participant)
    for b = 1:length(treatment)
        for c = 1:length(time)
            % load the header
            load([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.lw6'], '-mat') 

            % rename all events
            for e = 1:length(header.events)
                header.events(e).code = eventcode;
            end

            % save the header
            save([prefix ' ' participant{a} ' ' treatment{b} ' ' time{c} '.lw6'], 'header')
        end
    end
end
clear a b c e 

%% rename datasets
% parameters
participant = {'12' '13'};
treatment = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix_old = 'crop but fft-notchfilt prefilt prea dc reref ds art-sup ep dc EEG'; 
prefix_new = 'crop bandpass notch prefilt prea dc reref ds art-sup ep dc EEG'; 

% a = 1; b = 1; c = 1;
for a = 1:length(participant)
    for b = 1:length(treatment)
        for c = 1:length(time)
            % deal with the header
            load([prefix_old ' ' participant{a} ' ' treatment{b} ' ' time{c} '.lw6'], '-mat') 
            header.name = [prefix_new ' ' participant{a} ' ' treatment{b} ' ' time{c}];
            save([header.name '.lw6'], 'header')

            % deal with the data
            load([prefix_old ' ' participant{a} ' ' treatment{b} ' ' time{c} '.mat']) 
            save([header.name '.mat'], 'data')
        end
    end
end
clear a b c 