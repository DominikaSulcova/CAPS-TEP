%% session info
clear all; clc;

% parameters
timepoint = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'}; 
block = {'b1' 'b2'};

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
prefix = ['EEG ' session_info{1} ' ' session_info{2}];

% choose the folder with raw data
path = ['E:\Data\CAPS-TEP - data\Raw data\CAPS-TEP_' session_info{1}];
input_folder = uigetdir(path);

% specify EEG block subfolders in the imput folder
prompt = {'Admit following EEG blocks:'};
dlgtitle = 'EEG blocks';
dims = [1 35];
definput = {'[4:17]'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
eval(['folders = ' answer ';']); 
clear answer prompt dlgtitle dims definput             


%% import MEGA datasets
% import the datasets - blocks indicated by folders vector
counter = 1;
for a = 1:length(folders)
    % import the appropriate dataset
    [header, data] = EEG_import_MEGA(input_folder, folders(a));
                
    % create the name for the dataset
    if mod(a, 2) == 1
        dataset_name = [prefix ' ' timepoint{counter}  ' ' block{1}];
    else
        dataset_name = [prefix ' ' timepoint{counter}  ' ' block{2}];
        counter = counter + 1;
    end
    
    % create the first letswave history entry
    load('EEG_history_import.mat')
    EEG_history_import.configuration.parameters.input_folder  = input_folder;
    EEG_history_import.configuration.parameters.session_number  = a;   
    header.history(1) =  EEG_history_import;
    
    % save the data and the header as letswave files
    header.name = dataset_name;
    save([dataset_name '.mat'], 'data');
    save([dataset_name '.lw6'], 'header');
end

%% get rid of the first blind event 
counter = 1;
for b = 1:length(timepoint)*2 
    % choose the dataset name
    if mod(b, 2) == 1
        dataset_name = [prefix ' ' timepoint{counter}  ' ' block{1}];
    else
        dataset_name = [prefix ' ' timepoint{counter}  ' ' block{2}];
        counter = counter + 1;
    end
    
    % load the header 
    load([dataset_name '.lw6'],'-mat');  
    
    % ditch extra 'out' category
    T = struct2table(header.events);                        % creates a table of the structure 'events'
    T.code = categorical(T.code);                           % in order to work with strings --> convert to categorical
    sortedT = T(T.code == 'Stimulation', :); 
    sortedT.code = cellstr(sortedT.code);                   % turns the categorical data back to string cells
    
    % verify the number of events
    message_1 = [num2str(size(sortedT, 1)) ' events found in the file "' dataset_name '".'];
    disp(message_1)
    
    % choose the number of the blind event 
    prompt = {'Which event do you want to discard?'};
    dlgtitle = 'Blind event';
    dims = [1 50];
    definput = {'1'};
    blind_event = inputdlg(prompt,dlgtitle,dims,definput);
    blind_event = str2double(blind_event{1});
    clear prompt dlgtitle dims definput
    
    % discard the blind event 
    if  blind_event == 0                        % in case there is no blind event 
        message_2 = 'No event removed.';
        disp(message_2)
        header.events = table2struct(sortedT); 
        header.events = header.events';
    else
        message_2 = ['Removing blind event n.' num2str(blind_event) '.'];
        disp(message_2)
        T_crop = sortedT;                             
        T_crop(blind_event, :) = [];
        header.events = table2struct(T_crop); 
        header.events = header.events';
    end
    
    % save the header
    save([dataset_name '.lw6'], 'header');    
end 
