%% IMPORT EMG DATA TO LETSWAVE
% Written by Dominika for the CAPS-TEP project (2020)
% 
% 1) Loads raw EMG datasets saved in VHDR format
% 2) Keeps only event markers eith the target eventcode (for MEGA --> 's1')
%       - check for duplicate events (identical event latency) and discard
%       repeats
% 3) Saves in the current directory as .mat data + .lw6 header

%% session info
clear all; clc;

% parameters
timepoint = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
eventcode = 's1';

% enter dataset info
prompt = {'Subject number:'};
dlgtitle = 'Subject';
dims = [1 35];
definput = {'00'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% condition
answer = questdlg('Which solution is applied?', 'Experimental condition',...
    'capsaicin', 'vehicle', 'none');
switch answer
    case 'capsaicin'
        session_info{2} = 'caps';
    case 'vehicle'
        session_info{2} = 'ctrl';
end

% create a prefix
prefix = ['EMG ' session_info{1} ' ' session_info{2}]; 

% choose the folder with raw data
path = ['E:\Data\CAPS-TEP - data\Raw data\CAPS-TEP_' session_info{1} '\EMG'];
input_folder = uigetdir(path);
EMG_code = uigetfile([path '\*.vhdr'], 'Select file name for this session');
clear answer

%% import MEGA datasets
for a = 1:length(timepoint)
    % import the appropriate dataset
    filename = [input_folder '\' timepoint{a} '\' EMG_code];
    [header, data] = EMG_import_VHDR(filename);
                   
    % create the name for the dataset
    dataset_name = [prefix ' ' timepoint{a}];

    % create the first history entry
    load('EMG_history_import.mat')
    EMG_history_import.configuration.parameters.filenames =  filename;
    header.history(1) =  EMG_history_import;
    
    % keep only the events with target code 
    for b = 1:length(header.events)
        if strcmp(header.events(b).code, eventcode)
            index(b) = true;
        else
            index(b) = false;
        end
    end
    header.events = header.events(index);
    
    % remove repeated events (based on latency)
    events_unique = unique(extractfield(header.events, 'latency'));
    index_rep = [];
    for c = 1:length(events_unique)
        e = find(extractfield(header.events, 'latency') == events_unique(c));
        index_rep(end + 1) = e(1);
    end
    header.events = header.events(index_rep);
    
    % save the data and the header as letswave files
    header.name = dataset_name;
    save([dataset_name '.mat'], 'data');
    save([dataset_name '.lw6'], 'header');
end
clear a b c e events_unique index index_rep