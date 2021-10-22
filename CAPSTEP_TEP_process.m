%% CAPSTEP: PROCESS EeG DATA
% Written by Dominika for the CAPS-TEP project (2021)

% ----- Prepare data ----- 
% 1) load the data
%    - loads preprocessed data from letswave
%    - crops the data in predefined time window
%       - default [-0.05, 0.3]s
%   - saves to the global MATLAB file ('filename') --> 'CAPSTEP_data'
% 
% 2) export data for RAGU
%   - removes 'target' channel
%   - saves time-series in a .csv table, timepoint x channel (export folder)
%   - creates an .xyz file with the electrode montage

%% parameters
clear all; clc

% ----- adjustable parameters -----
% output
filename = 'CAPSTEP_prelim';

% dataset
participant = {'02' '03' '04' '05' '06' '07' '08' '10' '12' '13' '14' '16' '17'};
session = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
prefix = 'avg avgchan bl icfilt ica visual crop bandpass notch prefilt prea dc reref ds art-sup ep dc EEG'; 

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

% create output folders
folderpath = uigetdir;
if ~exist([folderpath '\CAPSTEP_TEP_export'])
    mkdir(folderpath, 'CAPSTEP_TEP_export')
end
if ~exist([folderpath '\CAPSTEP_TEP_figures'])
    mkdir(folderpath, 'CAPSTEP_TEP_figures')  
end

% header
load([prefix ' 02 caps baseline.lw6'], '-mat')
labels = {header.chanlocs.labels};

% visualization 
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

% colours
if ~exist([folderpath '\colours.mat'])
    for c = 1:length(session)
        colours(c, :) = uisetcolor(['Choose colour for ' session{c} ' data:']);
    end
    save([folderpath '\colours.mat'], 'colours')
    clear c
else
    load('colours.mat')
end

%% 1) extract individual data
% load data --> uncorrected ppTMS TEPs
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time) 
            % load individual dataset
            load([prefix ' ' participant{p} ' ' session{s} ' ' time{t} '.mat'])

            % append the data in the data matrix
            CAPSTEP_data(p, s, t, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                  
        end
    end
end
clear p s t data 
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist([folderpath '\' filename '.mat']) 
    save([folderpath '\' filename '.mat'], 'CAPSTEP_data');
else
    save([folderpath '\' filename '.mat'], 'CAPSTEP_data', '-append');
end

%% 2) export data for Ragu
% write text files for Ragu 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time) 
            % choose data to write, remove 'target' channel
            data = squeeze(CAPSTEP_data(p, s, t, 1:30, :))';

            % save as .csv               
            name = ['CAPSTEP_' participant{p} '_' session{s} '_' time{t} '.csv']; 
            writematrix(data, [folderpath '\CAPSTEP_TEP_export\' name])
        end
    end
end
clear p s t data name     

% create the montage file
name = [folderpath '\CAPSTEP_TEP_export\CAPSTEP_montage.xyz'];
fileID = fopen(name, 'a');
fprintf(fileID, '30\r\n');
for a = 1:30
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear name fileID a 















