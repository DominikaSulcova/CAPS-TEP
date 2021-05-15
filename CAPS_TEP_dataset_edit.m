%% rename categories
clear all, clc

% parameters
participant = {'08'};
session = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix = 'dc reref ds art-sup ep dc EEG';
event_name = 'Stimulation';


% rename loop 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)
            load([prefix ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'], '-mat')
            for e = 1:length(header.events)
                header.events(e).code = event_name;
            end
            save([prefix ' ' participant{p} ' ' session{s} ' ' time{t} '.lw6'], 'header')
        end
    end
end
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

