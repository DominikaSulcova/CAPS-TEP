%% IMPORT EMG DATA TO LETSWAVE
% 

%% session info
clear all; clc;

% parameters
timepoint = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};

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
% import the datasets - blocks indicated by folders vector
counter = 1;
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
      
    % save the data and the header as letswave files
    header.name = dataset_name;
    save([dataset_name '.mat'], 'data');
    save([dataset_name '.lw6'], 'header');
end