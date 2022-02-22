%% CAPSTEP: PROCESS EEG DATA
% Written by Dominika for the CAPS-TEP project (2021)
% 
% - each section of the script performs one data processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - when the manual input is needed (ICA, visual inspection), the process
%   is performed in the letswave 7 GUI and the corresponding section of the
%   script takes care of the logfile update and encodes the output
%   parameters in a structure 'CAPSTEP_parameters'
% 1) processing block 1
%     - re-reference to Cz
%     - remove mastoids - channels M1(29) and M2(30)
%     - re-reference to the common average
%     - segment into epochs [-1 2]s relative to the stimulus  
%     - remove DC + apply linear detrend
%     - replace the TMS artifact [-0.005 0.01]s - cubic interpolation
%     - downsample by 1/10 - final SR 2kHz
% 2) merge blocks 1 and 2
% 3) remove bad events
% 4) preliminary ICA matrix
%     - square matrix --> 30 components
% 5) preliminary ICA - encode results
% 6) processing block 2
%     - bandpass Butterworth filter [0.1 80]Hz
%     - notch FFT filter 50Hz
%     - signal cropped to [-0.75 0.75]s
% 7) csecond ICA matrix
%     - rectangular matrix --> number of components for each subject stored
%     in a vector 'components'
% 8) visual inspection - encode results


%% parameters
clear all
clc
% ------ participants -----
participant = [2:10, 12:14, 16:18, 20, 22];
% -------------------------

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
participant = [2:10, 12:14, 16:18, 20, 22];
condition = {'caps' 'ctrl'};
time = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
block = {'b1' 'b2'};
prefix = 'EEG';

% output
output_file = [folder_results '\CAPSTEP_parameters.mat' ];

%% 0) CAPSTEP parameters
% ----- section input -----
participant = [1:10, 12:14, 16:18, 20, 22];
% -------------------------
% check for an existing file
if ~exist(output_file) 
    CAPSTEP_parameters = struct;
else
    load(output_file);
end

% extract EMG identifiers
for p = 1:length(participant)
    % subject number
    CAPSTEP_parameters(participant(p)).subject_ID = participant(p);
    
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
        statement = ['CAPSTEP_parameters(participant(p)).EMG_ID.' condition{c} ' = EMG_ID;'];
        eval(statement)
    end
end
clear p c subject filename fileID nline EMG_ID statement

% save the output file
save(output_file, 'CAPSTEP_parameters');

%% 1) processing block 1
% ----- section input -----
suffix = {'reref' 'chan-select' 'reref' 'ep' 'dc' 'art-sup' 'ds'};
% participant = [];
% -------------------------
for p = 1:length(participant)
    for c = 1:length(condition)
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('block 1', filename)
    end
end
clear p c 

%% 2) merge block 1 and 2
% ----- section input -----
channel_n = 30;
% participant = [];
% -------------------------
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
           % load both datasets
           for b = 1:length(block)
               load([name ' ' char(time(t)) ' ' char(block(b)) '.mat']);
               load([name ' ' char(time(t)) ' ' char(block(b)) '.lw6'], '-mat');    
               statement = ['data_' num2str(b) ' = data; header_' num2str(b) ' = header;'];
               eval(statement);
           end

           % merge the data together
           data = cat(1, data_1, data_2); 

           % merge events and update in header
           T_1 = struct2table(header_1.events); 
           T_2 = struct2table(header_2.events); 
           T = [T_1; T_2];
           n = size(T, 1);
           T.epoch(1:end) = [1:n];
           header.events = table2struct(T);                  
           header.events = header.events';  

           % update the header info
           header.name = [name ' ' char(time(t))];
           if header_1.datasize(2) == channel_n & header_2.datasize(2) == channel_n
               header.datasize = [n, channel_n, 1, 1, 1, 6000];  
           else
               disp('ERROR: uneaven number of channels!')
               return
           end

           % save the merged dataset
           save([name ' ' char(time(t)) '.mat'], 'data');
           save([name ' ' char(time(t)) '.lw6'], 'header');
        end
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('merge', filename)
    end
end
clear p c t b name statement data_1 data_2 data...
    header_1 header_2 header T T_1 T_2 n channel_n

%% 3) remove bad events
% ----- section input -----
participant = 1;
% -------------------------
for p = 1:length(participant)
    for c = 1:length(condition)
        % ask for removed events
        prompt = {'Events removed:'};
        dlgtitle = 'Bad events';
        dims = [1 35];
        definput = {'[]'};
        answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
        eval(['removed = ' answer ';']); 
        clear answer prompt dlgtitle dims definput   
        
        % encode to the output structure        
        statement = ['CAPSTEP_parameters(participant(p)).bad_events.' condition{c} ' = removed;'];
        eval(statement)
        
        % identify events    
        if length(removed) > 0
            for r = 1:length(removed)
                statement = ['[removed_id(r, 1), removed_id(r, 2)] = find(CAPSTEP_parameters(1).EMG_ID.' condition{c} ' == removed(r));'];
                eval(statement)
            end
            values = unique(removed_id(:, 1));
            for u = 1:length(values)
                values_i = removed_id(removed_id(:, 1) == values(u), 2);
                removed_ID(u, 1) = {values(u)};
                removed_ID(u, 2) = {values_i};
            end
        else
            removed_ID = [];
        end
        clear r u values values_i removed_id
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('bad events', filename, 'removed', removed_ID)
    end
end
clear p c statement filename

% save the output file
save(output_file, 'CAPSTEP_parameters');

%% 4) preliminary ICA matrix
% ----- section input -----
suffix = 'prea';
prefix = 'ds art-sup dc ep reref chan-select reref EEG';
% participant = [];
% -------------------------
% calculate ICA matrix
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

%% 5) preliminary ICA 
% ----- section input -----
suffix = 'prefilt';
participant = 1;
% -------------------------
for p = 1:length(participant)
    for c = 1:length(condition)
        % ask for removed components
        prompt = {'Components removed:'};
        dlgtitle = 'Preliminary ICA';
        dims = [1 35];
        definput = {'[]'};
        answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
        eval(['components = ' answer ';']); 
        clear answer prompt dlgtitle dims definput  
        
        % encode to the output structure        
        statement = ['CAPSTEP_parameters(participant(p)).bad_events.' condition{c} ' = removed;'];
        eval(statement)
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('pre-ICA', filename, 'components', components)
        
        % rename
    end
end
clear p c filename

% update prefix
prefix = [suffix ' ' prefix];

% save the output file
save(output_file, 'CAPSTEP_parameters');

%% 6) processing block 2
% ----- section input -----
participant = [2:10, 12:14, 16:18, 20, 22];
suffix = {'bandpass' 'notch' 'crop'};
% -------------------------

%% 7) second ICA
% ----- section input -----
participant = [2:10, 12:14, 16:18, 20, 22];
components2run = 30 - components;
suffix = 'ica';
% -------------------------
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

        % compute the ICA matrix
        option = struct('ICA_mode', 2, 'algorithm', 1, 'suffix', suffix, 'is_save', 1);
        lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
        clear suffix
    end
end

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
            if length(removed) == 0
                fprintf(fileID, '4	no missed stimuli --> no faulty events discarded\r\n');
                fprintf(fileID, '\r\n');
                fclose(fileID);
            else
                fprintf(fileID, '4	faulty (missed) events were removed:\r\n');
                fprintf(fileID, '\r\n');
                fclose(fileID);
            end

        case 'pre-ICA'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'removed'));
                if ~isempty(a)
                    removed = varargin{a + 1};
                end
            end
            
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
    end
end
