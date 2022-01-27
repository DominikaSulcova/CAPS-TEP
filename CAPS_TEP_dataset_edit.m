%% rename categories
clear all, clc

% parameters
participant = {'03' '04' '05' '07' '08'};
session = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
block = {'b1' 'b2'};
event_name = 'Stimulation';

% loop through datasets
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)
            for b = 1:length(block)
                % load the header
                load(['EEG ' participant{p} ' ' session{s} ' ' time{t} ' ' block{b} '.lw6'], '-mat')
                
                % replace event code
                for e = 1:length(header.events)
                    header.events(e).code = event_name;
                end
                
                % save to letswave
                save(['EEG ' participant{p} ' ' session{s} ' ' time{t} ' ' block{b} '.lw6'], 'header')
            end
        end
    end
end
clear p s t 

%% rename datasets
clear all, clc

% parameters
participant = {'03' '04' '05' '06' '07' '08' '10' '14' '16' '17'};
session = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix_old = 'crop but fft-notchfilt prefilt prea dc reref ds art-sup ep dc EEG';
prefix_new = 'crop bandpass notch prefilt prea dc reref ds art-sup ep dc EEG';

% rename loop 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)
            load([prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.mat'])
            load([prefix_old ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'], '-mat')
            header.name = [prefix_new ' ' participant{p} ' ' session{s} ' ' time{t}];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header')
        end
    end
end

