%% CAPSTEP: PROCESS EEG DATA
% Written by Dominika for the CAPS-TEP project (2021)
% 
% - each section of the script performs one data processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - when the manual input is needed (ICA, visual inspection), the process
%   is performed in the letswave 7 GUI and the corresponding section of the
%   script takes care of the logfile update and encodes the output
%   parameters in a structure 'CAPSTEP_info'
% 
% 1) CAPSTEP info
%     - encode subject ID
%     - encode MEP IDs 
% 2) processing block 1
%     - re-reference to Cz
%     - remove mastoids - channels M1(29) and M2(30)
%     - re-reference to the common average
%     - segment into epochs [-1 2]s relative to the stimulus  
%     - remove DC + apply linear detrend
%     - replace the TMS artifact [-0.005 0.01]s - cubic interpolation
%     - downsample by 1/10 - final SR 2kHz
% 3) merge blocks 1 and 2
% 4) bad events
%     - asks for MEP IDs of bad events
%     - identifies their location, removes from the dataset
%     --> CAPSTEP info: encode removed events
% 5) preliminary ICA matrix
%     - square matrix --> 30 components
% 6) preliminary ICA
%     --> CAPSTEP info: encode removed ICs 
% 7) processing block 2
%     - bandpass Butterworth filter [0.1 80]Hz
%     - notch FFT filter 50Hz
%     - signal cropped to [-0.75 0.75]s
% 8) second ICA matrix
%     - rectangular matrix --> number of components for each subject stored
%     in a vector 'components'
% 9) second ICA
%     --> CAPSTEP info: encode removed ICs  
% 10) visual inspection
%     --> CAPSTEP info: encode removed epochs 
% 11) target channel
%     - create a new channel by averaging 3 closest electrodes
% 12) processing block 3
%     - baseline correction: subtraction of baseline segment [-0.2 -0.005]s
%     - average across epochs


%% parameters
clear all
clc

% choose the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% choose the folder with processed data 
path = 'E:\Data\CAPS-TEP - data\Processed data';
folder_output = uigetdir(path, 'Coose the output folder');

% choose the results folder  
path = 'E:\UCL\O365G-NOCIONS - CAPS-TEP - CAPS-TEP\Results';
folder_results = uigetdir(path, 'Coose the results folder');

% choose the folder with logfiles
path = 'E:\UCL\O365G-NOCIONS - CAPS-TEP - CAPS-TEP\Logfiles';
folder_logfiles = uigetdir(path, 'Coose the logfile folder');
clear path 

% dataset
condition = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
block = {'b1' 'b2'};
prefix = 'EEG';

% output
output_file = [folder_results '\CAPSTEP_info.mat' ];

%% 1) CAPSTEP info: encode MEP IDs 
% ----- section input -----
% participant = [];
% -------------------------
% check for an existing file
if ~exist(output_file) 
    CAPSTEP_info = struct;
else
    load(output_file);
end

% loop through datasets
for p = 1:length(participant)
    % subject number
    CAPSTEP_info(participant(p)).subject_ID = participant(p);
    
    % identify EMG events
    for c = 1:length(condition)
        % determine the filename
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        
        % get the line number
        fileID = fopen(filename, 'r');
        nline = getline(fileID, '	- extracted from Visor2');
        
        % extract the content
        fileID = fopen(filename, 'r');
        EMG_ID = get_emgid(fileID, nline);

        % encode to the output structure        
        statement = ['CAPSTEP_info(participant(p)).EMG_ID.' condition{c} ' = EMG_ID;'];
        eval(statement)
    end
end
clear p c subject filename fileID nline EMG_ID statement

% save the output file
save(output_file, 'CAPSTEP_info');

%% 2) processing block 1
% ----- section input -----
suffix = {'reref' 'chan-select' 'reref' 'ep' 'dc' 'art-sup' 'ds'};
electrodes = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2','CP1','CP2','FC5','FC6','CP5','CP6','M1','M2','C1','C2'};
% participant = [];
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        for t = 1:length(time)
            for b = 1:length(block)
                % get the dataset
                disp([subject ' ' condition{c} ' ' time{t} ' ' block{b}])
                dataset = [folder_output '\' prefix ' ' subject ' ' condition{c} ' ' time{t} ' ' block{b} '.lw6'];
                option = struct('filename', dataset);
                lwdata = FLW_load.get_lwdata(option);

                % re-reference to Cz
                disp('Re-referencing to Cz...')
                option = struct('reference_list', {{'Cz'}}, 'apply_list', {electrodes}, 'suffix', suffix{1}, 'is_save', 0);
                lwdata = FLW_rereference.get_lwdata(lwdata, option);
                
                % remove mastoids
                disp('Removing mastoids...')
                option = struct('type', 'channel', 'items', {electrodes([1:28, 31, 32])}, 'suffix', suffix{2}, 'is_save', 0);
                lwdata = FLW_selection.get_lwdata(lwdata, option);
                
                % re-reference to common average
                disp('Re-referencing to common average...')
                option = struct('reference_list', {electrodes([1:28, 31, 32])}, 'apply_list', {electrodes([1:28, 31, 32])},...
                    'suffix', suffix{3}, 'is_save', 0);
                lwdata = FLW_rereference.get_lwdata(lwdata, option);
                
                % segment
                disp('Epoching from -1 to 2 s relative to the event...')
                option = struct('event_labels', {{'Stimulation'}}, 'x_start', -1, 'x_end', 2, 'x_duration', 3, ...
                    'suffix', suffix{4}, 'is_save', 0);
                lwdata = FLW_segmentation.get_lwdata(lwdata, option);
                
                % remove DC + linear detrend
                disp('Removing DC and applying linear detrend...')
                option = struct('linear_detrend', 1, 'suffix', suffix{5}, 'is_save', 0);
                lwdata = FLW_dc_removal.get_lwdata(lwdata,option);
                
                % artifact removal, interpolation
                disp('Interpolating the TMS artifact from -0.005 to 0.01 s...')
                [lwdata.header, lwdata.data, message_string] = RLW_suppress_artifact_event(lwdata.header, lwdata.data,...
                    'xstart', -0.005, 'xend', 0.01, 'event_code', 'Stimulation', 'interp_method', 'pchip');
                                
                % downsample
                disp('Downsampling...')
                option = struct('x_dsratio', 10,'suffix', [suffix{7} ' ' suffix{6}], 'is_save', 1);
                lwdata = FLW_downsample.get_lwdata(lwdata, option);
            end
        end
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('block 1', filename)
    end
end
clear p c subject t b dataset option lwdata

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s

%% 3) merge block 1 and 2
% ----- section input -----
channel_n = 30;
% participant = [];
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        % loop through timepoints
        for t = 1:length(time)
            % load the dataset
            for b = 1:length(block)
                dataset{b} = [folder_output '\' prefix ' ' subject ' ' condition{c} ' ' time{t} ' ' block{b} '.lw6'];
            end
            option = struct('filename', {dataset});
            lwdataset = FLW_load.get_lwdataset(option);
        
            % merge epochs 
            option = struct('type', 'epoch', 'suffix', '', 'is_save', 0);
            lwdata = FLW_merge.get_lwdata(lwdataset, option);
            
            % rename, save
            lwdata.header.name = [prefix ' ' subject ' ' condition{c} ' ' time{t}];
            CLW_save(lwdata)       
        end
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('merge', filename)
    end
end
clear p c t b subject dataset option lwdataset lwdata filename

%% 4) remove bad events
% ----- section input -----
% participant = [];
% -------------------------
% load output file
load(output_file);

% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
                
        % ask for removed events
        prompt = {'Events to remove:'};
        dlgtitle = sprintf('Subject %d, %s: bad events', participant(p), condition{c});
        dims = [1 60];
        definput = {'[]'};
        answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
        eval(['events2remove = ' answer ';']); 
        clear answer prompt dlgtitle dims definput   
        
        % identify events   
        removed_ID = {};
        if length(events2remove) > 0
            % identify MEP position
            for e = 1:length(events2remove)
                statement = ['[removed_id(e, 1), removed_id(e, 2)] = find(CAPSTEP_info(participant(p)).EMG_ID.' condition{c} ' == events2remove(e));'];
                eval(statement)
            end
            
            % transform into a cell array with one row per timepoint
            values = unique(removed_id(:, 1));
            for v = 1:length(values)
                values_i = removed_id(removed_id(:, 1) == values(v), 2)';
                removed_ID(v, 1) = {values(v)};
                removed_ID(v, 2) = {values_i};
            end
        end
        clear e v values values_i removed_id 
        
        % remove the events from appropriate datasets
        kept_ID = {};
        for d = 1:size(removed_ID, 1)
            % load the dataset
            dataset = [folder_output '\' prefix ' ' subject ' ' condition{c} ' ' time{removed_ID{d, 1}} '.lw6'];
            option = struct('filename', dataset);
            lwdata = FLW_load.get_lwdata(option);
            
            % remove events
            kept_id = [1:80];
            kept_id(removed_ID{d, 2}) = [];
            for k = 1:length(kept_id)
                kept_ID{k} = num2str(kept_id(k));
            end
            option = struct('type', 'epoch', 'items', {kept_ID}, 'suffix', '', 'is_save', 1);
            lwdata = FLW_selection.get_lwdata(lwdata, option);
        end
        clear d k kept_id 
        
        % encode to the output structure        
        statement = ['CAPSTEP_info(participant(p)).bad_events.' condition{c} ' = events2remove;'];
        eval(statement)
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('bad events', filename, 'removed', removed_ID)
    end
end
clear p c subject statement filename dataset option lwdata events2remove kept_ID removed_ID

% save the output file
save(output_file, 'CAPSTEP_info');

%% 5) preliminary ICA matrix
% ----- section input -----
suffix = 'prea';
% participant = [];
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        % load the dataset
        for t = 1:length(time)
            dataset{t} = [folder_output '\' prefix ' ' subject ' ' condition{c} ' ' time{t} '.lw6'];
        end
        option = struct('filename', {dataset});
        lwdataset = FLW_load.get_lwdataset(option);

        % compute the ICA matrix
        option = struct('ICA_mode', 1, 'algorithm', 1, 'suffix', suffix, 'is_save', 1);
        lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    end
end
clear p c t subject dataset option lwdataset

% update prefix
prefix = [suffix ' ' prefix];

%% 6) preliminary ICA 
% ----- section input -----
suffix = 'prefilt';
participant = [1:10, 12:14, 16:18, 20, 22:24];
% -------------------------
% load output file
load(output_file);

% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        % load one of the datasets
        dataset = [folder_output '\sp_filter ' prefix ' ' subject ' ' condition{c} ' ' time{1} '.lw6'];
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);
            
        % extract numbers of removed components
        components = lwdata.header.history(2).option.remove_idx;  
        
        % encode to the output structure        
        statement = ['CAPSTEP_info(participant(p)).preICA.' condition{c} ' = components;'];
        eval(statement)
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('preICA', filename, 'components', components)
        
        % rename, save
        for t = 1:length(time)
            % load the dataset
            dataset = [folder_output '\sp_filter ' prefix ' ' subject ' ' condition{c} ' ' time{1} '.lw6'];
            option = struct('filename', dataset);
            lwdata = FLW_load.get_lwdata(option);
            
            % rename
            lwdata.header.name = [suffix ' ' prefix ' ' subject ' ' condition{c} ' ' time{t}];
            CLW_save(lwdata)    
        end       
    end
end
clear p c subject filename dataset option lwdata t

% update prefix
prefix = [suffix ' ' prefix];

% save the output file
save(output_file, 'CAPSTEP_info');

%% 7) processing block 2
% ----- section input -----
suffix = {'bandpass' 'notch' 'crop'};
% participant = [];
% -------------------------
% loop through datasets

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s

%% 8) second ICA matrix
% ----- section input -----
suffix = 'ica';
% participant = [];
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        % load the dataset
        for t = 1:length(time)
            dataset{t} = [folder_output '\' prefix ' ' subject ' ' condition{c} ' ' time{t} '.lw6'];
        end
        option = struct('filename', {dataset});
        lwdataset = FLW_load.get_lwdataset(option);
        clear t subject dataset
        
        % determine number of components to run
        components = numel(lwdataset(1).header.history(2).option.remove_idx);
        components2run = 30 - components;

        % compute the ICA matrix
        option = struct('ICA_mode', 2, 'algorithm', 1, 'suffix', suffix, 'is_save', 1);
        lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    end
end

% update prefix
prefix = [suffix ' ' prefix];

%% 9) second ICA
% ----- section input -----
suffix = 'ica_filt';
% participant = [];
% -------------------------
% load output file
load(output_file);

% loop through datasets

% update prefix
prefix = [suffix ' ' prefix];

% save the output file
save(output_file, 'CAPSTEP_info');

%% 10) visual inspection
% ----- section input -----
suffix = 'visual';
% participant = [];
% -------------------------
% load output file
load(output_file);

% loop through datasets

% update prefix
prefix = [suffix ' ' prefix];

% save the output file
save(output_file, 'CAPSTEP_info');

%% 11) target channel
% ----- section input -----
% participant = [];
% -------------------------
% load output file
load(output_file);

% loop through datasets

% save the output file
save(output_file, 'CAPSTEP_info');

%% 12) processing block 3
% ----- section input -----
suffix = {'bl' 'avg'};
% participant = [];
% -------------------------
% loop through datasets

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s

%% functions
function nline = getline(fileID, indice)
    tline = fgetl(fileID);
    nline = 1;
    while ~strcmp(tline, indice)
        tline = fgetl(fileID);
        nline = nline + 1;
    end
end
function EMG_ID = get_emgid(fileID, nline)
    pos = [4, 6, 10, 12, 16, 18, 22, 24, 28, 30, 34, 36, 40, 42];
    textfile = textscan(fileID, '%s', 'Headerlines', nline);
    for a = 1:length(pos)
        n_read = cell2mat(textfile{1, 1}(pos(a)));  
        delim = find(n_read == '-');
        n_start(a) = str2double(n_read(1:delim-1));
        n_end(a) = str2double(n_read(delim+1:end));
    end
    EMG_ID = [];
    for b = 1:2:length(n_start)
        values = [n_start(b):n_end(b), n_start(b+1):n_end(b+1)];
        EMG_ID(end+1, 1:length(values)) = values;
    end
end
function logfile_entry(entry, filename, varargin)
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
    switch entry
        case 'block 1'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, 'DATA CLEANING\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, 'TEPs over M1\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, '1	dataset preparation: used a matlab script ''CAPSTEP_EEG_import.m''\r\n');
            fprintf(fileID, '		- data automatically imported\r\n');
            fprintf(fileID, '		- ''Out'' category discarded, ''Stimulation'' category kept\r\n');
            fprintf(fileID, '		- checked for missing/extra events, first blind event discarded\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '2	first step preprocessing: used a matlab script ''CAPSTEP_EEG_process.m''\r\n');
            fprintf(fileID, '		- electrodes were assigned standardized 10/20 labels\r\n');
            fprintf(fileID, '		- rereferenced to Cz\r\n');
            fprintf(fileID, '		- mastoids removed\r\n');
            fprintf(fileID, '		- rereferenced to the common average\r\n');
            fprintf(fileID, '		- segmentation relative to events --> epochs [-1 2]s, epoch size 60000 bins\r\n');
            fprintf(fileID, '		- DC removal + linear detrend\r\n');
            fprintf(fileID, '		- TMS artifact suppression --> [-0.005 0.01]s\r\n');
            fprintf(fileID, '		- downsample by 1/10\r\n');
            fprintf(fileID, '	--> file prefix: ds art-sup dc ep reref chan-select reref\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);

        case 'merge'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '3	blocks 1 and 2 were merged together and renamed: used a matlab script ''CAPSTEP_EEG_import.m''\r\n');
            fprintf(fileID, '	--> ds art-sup dc ep reref chan-select reref %s %s timepoint\r\n', filename(end-10 : end-9), filename(end-7 : end-4)); 
            fprintf(fileID, '\r\n');
            fclose(fileID);

        case 'bad events'
            a = find(strcmpi(varargin, 'removed'));
            if ~isempty(a)
                removed = varargin{a + 1};
            end
            
            fileID = fopen(filename, 'a');
            if size(removed) == 0
                fprintf(fileID, '4	no missed stimuli --> no faulty events discarded\r\n');
                fprintf(fileID, '\r\n');
                fclose(fileID);
            else
                fprintf(fileID, '4	faulty (missed) events were removed:\r\n');
                for r = 1:size(removed, 1)
                    fprintf(fileID, '		- n. %s in %s\r\n', num2str(removed{r, 2}), time{removed{r, 1}});
                end
                fprintf(fileID, '\r\n');
                fclose(fileID);
            end

        case 'preICA'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'components'));
                if ~isempty(a)
                    components = varargin{a + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '5	preliminary ICA was computed to get rid of the major decay artifact\r\n');
            fprintf(fileID, '		- squared matrix --> 30 ICs\r\n');
            fprintf(fileID, '		- removed component(s): IC %s\r\n', num2str(components));
            fprintf(fileID, '	--> file prefix: prefilt prea');
            fprintf(fileID, '\r\n');
            fclose(fileID);
    end
end
