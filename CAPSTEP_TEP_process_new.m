%% CAPS-TEP: TEP data analysis
% ------------------------------------------------------------------------
% author:   Dominik Sulcova
%           MSH Hamburg, Germany
% created:  January 2025
% ------------------------------------------------------------------------
% data:     raw EEG recordings (NeurOne, Bittium)
%           - SR 20 kHz
%           - referenced to mastoids
%           - ground at AFz
% script:
%   1) data import & first pre-processing              
%       - loads raw data from the subject folder --> structure 'dataset'
% 
%% parameters
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Choose the input folder');              % raw data --> hard-drive or local folder
folder.processed = uigetdir(pwd, 'Choose the data folder');         % processed data --> local folder 
folder.output = uigetdir(pwd, 'Choose the output folder');          % output folder --> OneDrive folder
cd(folder.processed)

% output
study = 'CAPSTEP';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);

% load output structure
fprintf('loading the info structure...\n')
if exist(output_file) == 2
    output_vars = who('-file', output_file);
    if ismember('TEP_new', output_vars)
        load(output_file, 'TEP_new')
    else
        TEP_new = struct;
        save(output_file, 'TEP_new','-append')
    end
else
    TEP_new = struct; 
    save(output_file, 'TEP_new')
end

% current participant
prompt = {'current subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});

% visualization
figure_counter = 1;
clear output_vars prompt dlgtitle dims definput input

%% 1) data import & first pre-processing
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
params.suffix = {'ep' 'ds' 'dc' 'art_interp'};
params.eventcode = {'TMS'};
params.epoch = [-1.5 1.5];
params.downsample = 20;
params.artifact_interp = [-0.005 0.01];
params.artifact_method = 'pchip';
% ------------------------- 
fprintf('section 1: load & pre-process EEG data\n')

% cycle through sessions
for a = 1:length(params.condition)
    % add letswave 6 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 6']));

    % provide update
    fprintf('subject %d - %s session:\n', subject_idx, params.condition{a})

    % identify dataset
    dataset_idx = strcmp({TEP_new(subject_idx).recording.dataset.session}, 'pain');
    files2import = TEP_new(subject_idx).recording.dataset(dataset_idx);

    % identify import folder 
    for b = 1:length(TEP_new(subject_idx).session)
        if strcmp(TEP_new(subject_idx).session(b).condition, params.condition{a})
            params.folder = sprintf('%s\\%s\\%s', folder.raw, TEP_new(subject_idx).ID, ...
                TEP_new(subject_idx).session(b).date);
        end
    end

    % check avaiable datasets
    data2import = dir(params.folder);
    if ~isempty(data2import)
        data_idx = logical([]);
        for c = 1:length(data2import)
            if isempty(str2num(data2import(c).name))
                data_idx(c) = true;
            else
                data2import(c).block = str2num(data2import(c).name);
                data_idx(c) = false;
            end
        end
        data2import(data_idx) = [];
        [~, file_idx] = sort([data2import.block]);
        data2import = data2import(file_idx);
        fprintf('%d datasets found in the directory.\n', length(data2import))
    else
        error('ERROR: No datasets found in the directory.\n')
    end

    % check number of datasets
    if length(data2import) ~= length(files2import)
        error('ERROR: This does not match with expected number of datasets (%d)\n.Please verify manually.\n', length(data2import))
    end

    % load datasets
    dataset(a).condition = params.condition{a};
    fprintf('Loading:\n')
    for d = 1:length(data2import)
        % provide update
        fprintf('%s:', files2import(d).name)
    
        % import the dataset
        [dataset(a).raw(d).header, dataset(a).raw(d).data, ~] = RLW_import_MEGA(data2import(d).folder, data2import(d).block);

        % rename in the header
        dataset(a).raw(d).header.name = files2import(d).name;
    end  
    fprintf('\n')

    % add letswave 7 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 7']));

    % segment and downsample 
    fprintf('Pre-processing:\n')
    for d = 1:length(dataset(a).raw)
        % provide update
        fprintf('%s:\n', files2import(d).name)

        % select data
        lwdata.header = dataset(a).raw(d).header;
        lwdata.data = dataset(a).raw(d).data; 

        % assign electrode coordinates
        fprintf('assigning electrode coordinates... ')
        option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
            'suffix', '', 'is_save', 0);
        lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(1).process = sprintf('electrode coordinates assigned');
            TEP_new(subject_idx).processing(1).params.layout = sprintf('standard 10-20-cap81');
            TEP_new(subject_idx).processing(1).date = sprintf('%s', date);
        end

        % identify eventcode
        fprintf('checking events... ') 
        if d == 1
            params.eventcode_old = unique({lwdata.header.events.code});
            prompt = 'Which eventcode corresponds to the TMS stimulus?';
            definput = '';
            for e = 1:length(params.eventcode_old)
                if e == length(params.eventcode_old)
                    definput = [definput params.eventcode_old{e}];
                else
                    definput = [definput params.eventcode_old{e} '\' ];
                end
            end
            dlgtitle = sprintf('subject %d - %s session: eventcodes', subject_idx, params.condition{a});
            dims = [1 40];
            params.eventcode_old = inputdlg(prompt,dlgtitle,dims,{definput});
            clear prompt dlgtitle dims definput 
        else
            if ~ismember(params.eventcode_old, unique({lwdata.header.events.code}))
                params.eventcode_old1 = params.eventcode_old;
                params.eventcode_old = unique({lwdata.header.events.code});
                prompt = sprintf(['Eventcode ''%s'' was not found in this dataset.\n' ...
                    'Which eventcode corresponds to the TMS stimulus?'], params.eventcode_old1{1});
                definput = '';
                for e = 1:length(params.eventcode_old)
                    if e == length(params.eventcode_old)
                        definput = [definput params.eventcode_old{e}];
                    else
                        definput = [definput params.eventcode_old{e} '\' ];
                    end
                end
                dlgtitle = sprintf('%s',  files2import(d).name);
                dims = [1 40];
                params.eventcode_old = inputdlg(prompt,dlgtitle,dims,{definput});
                clear prompt dlgtitle dims definput 
            end
        end

        % re-label and filter events
        event_idx = logical([]);
        for e = 1:length(lwdata.header.events)
            if strcmp(lwdata.header.events(e).code, params.eventcode_old{1})
                lwdata.header.events(e).code = params.eventcode{1};
                event_idx(e) = false; 
            else
                event_idx(e) = true; 
            end
        end
        lwdata.header.events(event_idx) = [];

        % remove first stimulus if the latency corresponds
        if lwdata.header.events(2).latency - lwdata.header.events(1).latency > 9.8
            lwdata.header.events(1) = [];    
        end
        fprintf('%d TMS stimuli found.\n', length(lwdata.header.events))

        % segment
        fprintf('segmenting ... ')
        option = struct('event_labels', params.eventcode, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            TEP_new(subject_idx).processing(2).process = sprintf('segmented to pre-stimulus epochs');
            TEP_new(subject_idx).processing(2).params.limits = params.epoch;
            TEP_new(subject_idx).processing(2).suffix = params.suffix{1};
            TEP_new(subject_idx).processing(2).date = sprintf('%s', date);
        end
        TEP_new(subject_idx).processing(2).params.epochs((a-1)*length(params.timepoint) + d) = size(lwdata.data, 1);

        % downsample 
        fprintf('downsampling... ')
        option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_downsample.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            TEP_new(subject_idx).processing(3).process = sprintf('downsampled');
            TEP_new(subject_idx).processing(3).params.ratio = params.downsample;
            TEP_new(subject_idx).processing(3).params.fs_orig = 1/lwdata.header.xstep * params.downsample;
            TEP_new(subject_idx).processing(3).params.fs_final = 1/lwdata.header.xstep;
            TEP_new(subject_idx).processing(3).suffix = params.suffix{2};
            TEP_new(subject_idx).processing(3).date = sprintf('%s', date);
        end

        % remove DC + linear detrend
        fprintf('removing DC and applying linear detrend.\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            TEP_new(subject_idx).processing(4).process = sprintf('DC + linear detrend on ERP epochs');
            TEP_new(subject_idx).processing(4).suffix = params.suffix{3};
            TEP_new(subject_idx).processing(4).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).raw(d).header = lwdata.header;
        dataset(a).raw(d).data = lwdata.data; 
    end
    fprintf('Done.\n\n')

    % add letswave 6 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 6']));

    % concatenate timepoints, interpolate TMS artifact and save for letswave
    fprintf('concatenating blocks, interpolating TMS artifact ... ') 
    for b = 1:length(params.timepoint)
        % identify blocks, check their number
        concat_idx = [];
        for d = 1:length(dataset(a).raw)
            if contains(dataset(a).raw(d).header.name, params.timepoint{b}) 
                concat_idx(end + 1) = d;
            end
        end
        if length(concat_idx) > 2
        elseif length(concat_idx) == 1
            fprintf('ATTENTION: Only 1 block was found for the %s dataset.\n', params.timepoint{b})
        elseif length(concat_idx) == 0
            fprintf('ATTENTION: No %s dataset was found!\n', params.timepoint{b})
            continue
        end

        % merge epochs of all identified datasets
        header = dataset(a).raw(concat_idx(1)).header;
        data = dataset(a).raw(concat_idx(1)).data;
        header.name = header.name(1:end-3); 
        if length(concat_idx) == 2
            for e = 1:size(dataset(a).raw(concat_idx(2)).data, 1)
                data(end+1, :, :, :, :, :) = dataset(a).raw(concat_idx(2)).data(e, :, :, :, :, :);
            end
        end
        header.datasize = size(data);

        % interpolate TMS artifact
        [header, data, ~] = RLW_suppress_artifact_event(header, data, ...
            'xstart', params.artifact_interp(1), 'xend', params.artifact_interp(2), ...
            'event_code', params.eventcode, 'interp_method', params.artifact_method); 
        header.name = [params.suffix{4} ' ' header.name];
        if d == 1 && c == 1
            INFO(subject_idx).processing(5).process = sprintf('TMS artifact interpolated');
            INFO(subject_idx).processing(5).params.limits = params.artifact_interp;
            INFO(subject_idx).processing(5).params.method = params.artifact_method;
            INFO(subject_idx).processing(5).suffix = params.suffix{4};
            INFO(subject_idx).processing(5).date = sprintf('%s', date);
        end        

        % save for letswave
        save([header.name '.lw6'], 'header')
        save([header.name '.mat'], 'data')

        % update dataset
        dataset(a).processed(b).header = header;
        dataset(a).processed(b).data = data;
    end
    fprintf('done.\n')
    fprintf('\n')
end

% save and continue
save(output_file, 'TEP_new','-append')
clear params a b c d e data2import dataset_idx files2import data_idx file_idx filename ...
    option lwdata event_idx epoch_idx concat_idx data header
fprintf('section 1 finished.\n')
