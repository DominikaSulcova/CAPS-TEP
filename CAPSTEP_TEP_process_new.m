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
%   1) ...
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

% dataset
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};

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

%% 1) import & pre-process continuous data
% ----- section input -----
params.suffix = {'crop' 'ds' 'dc'};
params.eventcode = 'TMS';
params.crop_margin = 5;
params.downsample = 20;
% ------------------------- 
fprintf('section 1: import & pre-process continuous data\n')

% load and pre-process continuous data
for a = 1:length(params.condition)
    % add letswave 6 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 6']));

    % provide update
    fprintf('subject %d - %s session:\n', subject_idx, params.condition{a})

    % identify dataset
    dataset_idx = strcmp({TEP_new(subject_idx).recording.dataset.session}, params.condition{a});
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
    fprintf('loading:\n')
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

    % pre-process 
    fprintf('pre-processing:\n')
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
                lwdata.header.events(e).code = params.eventcode;
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

        % crop 
        fprintf('cropping... ')
        params.crop(1) = lwdata.header.events(1).latency - params.crop_margin;
        params.crop(2) = lwdata.header.events(end).latency + params.crop_margin;
        option = struct('xcrop_chk', 1, 'xstart', params.crop(1), 'xend', params.crop(2), ...
            'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            TEP_new(subject_idx).processing(2).process = 'continuous data cropped';
            TEP_new(subject_idx).processing(2).params.start = params.crop(1);
            TEP_new(subject_idx).processing(2).params.end = params.crop(2);
            TEP_new(subject_idx).processing(2).params.margin = params.crop_margin;
            TEP_new(subject_idx).processing(2).suffix = params.suffix{1};
            TEP_new(subject_idx).processing(2).date = sprintf('%s', date);
        end

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
        fprintf('DC + linear detrend...\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 1);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            TEP_new(subject_idx).processing(4).process = sprintf('DC + linear detrend on continuous data');
            TEP_new(subject_idx).processing(4).suffix = params.suffix{3};
            TEP_new(subject_idx).processing(4).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).raw(d).header = lwdata.header;
        dataset(a).raw(d).data = lwdata.data; 
    end
    fprintf('done.\n\n')
end

% save and continue
save(output_file, 'TEP_new','-append')
clear a b c d e f data2import dataset_idx files2import data_idx data2load file_idx filename ...
    option lwdata event_idx epoch_idx concat_idx data header open fig_all
fprintf('section 1 finished.\n')

%% 2) segment and pre-process for visual inspection
% ----- section input -----
params.prefix = 'dc ds crop';
params.suffix = {'no_mastoid' 'ep' 'dc' 'art_interp' 'preprocessed'};
params.epoch = [-1.5 1.5];
params.artifact_interp = [-0.005 0.01];
params.artifact_method = 'pchip';
% ------------------------- 
fprintf('section 2: pre-processing for visual inspection\n')

% load output structure if needed 
if exist('TEP_new') ~= 1
    load(output_file, 'TEP_new')
end

% load dataset if needed
fprintf('loading dataset...\n')
if exist('dataset') ~= 1
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    if length(data2load) == length(params.condition) * length(params.timepoint) * length(params.block) * 2
        dataset = reload_dataset(data2load, params.condition, 'raw');
    else
        error(sprintf('ERROR: Wrong number of datasets (%d) found in the directory!', length(data2load)/2))
    end
end

% segment and pre-process
for a = 1:length(params.condition)
    for d = 1:length(dataset(a).raw)
        % provide update
        fprintf('%s:\n', dataset(a).raw(d).header.name)

        % add letswave 7 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 7']));

        % select data
        lwdata.header = dataset(a).raw(d).header;
        lwdata.data = dataset(a).raw(d).data; 

        % remove mastoids
        fprintf('removing mastoids ... ')
        channel_all = {lwdata.header.chanlocs.labels};
        channel_mask = cellfun(@(x) strcmp(x, 'M1') || strcmp(x, 'M2'), channel_all);
        channels2keep = channel_all(~channel_mask);
        option = struct('type', 'channel', 'items', {channels2keep}, 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_selection.get_lwdata(lwdata, option);
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(5).process = 'mastoid channels removed';
            TEP_new(subject_idx).processing(5).params.channels_kept = channels2keep;
            TEP_new(subject_idx).processing(5).suffix = params.suffix{1};
            TEP_new(subject_idx).processing(5).date = sprintf('%s', date);
        end

        % segment
        fprintf('segmenting ... ')
        option = struct('event_labels', params.eventcode, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(6).process = sprintf('segmented to pre-stimulus epochs');
            TEP_new(subject_idx).processing(6).params.limits = params.epoch;
            TEP_new(subject_idx).processing(6).suffix = params.suffix{2};
            TEP_new(subject_idx).processing(6).date = sprintf('%s', date);
        end

        % remove DC + linear detrend
        fprintf('DC + linear detrend ... ')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(7).process = sprintf('DC + linear detrend on epoched data');
            TEP_new(subject_idx).processing(7).suffix = params.suffix{3};
            TEP_new(subject_idx).processing(7).date = sprintf('%s', date);
        end

        % add letswave 6 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 6']));
        
        % interpolate TMS artifact
        fprintf('interpolating TMS artifact...\n')
        header = lwdata.header;
        data = lwdata.data;
        [header, data, ~] = RLW_suppress_artifact_event(header, data, ...
            'xstart', params.artifact_interp(1), 'xend', params.artifact_interp(2), ...
            'event_code', params.eventcode, 'interp_method', params.artifact_method); 
        header.name = [params.suffix{4} ' ' header.name];
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(8).process = sprintf('TMS artifact interpolated');
            TEP_new(subject_idx).processing(8).params.limits = params.artifact_interp;
            TEP_new(subject_idx).processing(8).params.method = params.artifact_method;
            TEP_new(subject_idx).processing(8).suffix = params.suffix{4};
            TEP_new(subject_idx).processing(8).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).raw(d).header = header;
        dataset(a).raw(d).data = data; 
    end
end
fprintf('done.\n')
fprintf('\n')

% concatenate blocks into timepoints
fprintf('concatenating blocks and saving ... ') 
for a = 1:length(params.condition)
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
        if length(concat_idx) == 2
            for e = 1:size(dataset(a).raw(concat_idx(2)).data, 1)
                data(end+1, :, :, :, :, :) = dataset(a).raw(concat_idx(2)).data(e, :, :, :, :, :);
            end
        end    
    
        % update header
        header.datasize = size(data);
        header.name = sprintf('%s %s %s %s %s', params.suffix{5}, study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b});
    
        % update dataset
        dataset(a).processed(b).header = header;
        dataset(a).processed(b).data = data;
    
        % save for letswave        
        save([header.name '.lw6'], 'header')
        save([header.name '.mat'], 'data')
    end
end
fprintf('done.\n')
fprintf('\n')

% open letswave for visual check
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    letswave
end

% save and continue
save(output_file, 'TEP_new','-append')
clear data2load a b d f lwdata option channel_all channel_mask channels2keep header data concat_idx fig_all open
fprintf('section 2 finished.\n')

%% 3) re-reference and export to EEGLAB
% ----- section input -----
params.prefix = 'preprocessed';
params.interp_chans = 6;
params.ref = 'averef';
% -------------------------
fprintf('section 3: exporting to EEGLAB\n')

% load output structure if needed 
if exist('TEP_new') ~= 1
    fprintf('loading output structure...\n')
    load(output_file, 'TEP_new')
end

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    if length(data2load) == length(params.condition) * length(params.timepoint) * 2
        dataset = reload_dataset(data2load, params.condition, 'processed');
    else
        error(sprintf('ERROR: Wrong number of datasets (%d) found in the directory!', length(data2load)/2))
    end
end

% interpolate channels if needed
params.labels = {dataset(1).processed(1).header.chanlocs.labels};
for c = 1:length(params.condition)
    prompt{c} = sprintf('%s session:', params.condition{c});
    definput{c} = strjoin(params.labels, ' ');
end
dlgtitle = 'channel interpolation';
dims = [1 120];
answer = inputdlg(prompt,dlgtitle,dims,definput);
for a = 1:length(answer)
    if ~isempty(answer{a})
        % identify channels to interpolate
        chans2interpolate = split(answer{a}, ' ');
    
        % interpolate identified channels in all datsets
        for c = 1:length(chans2interpolate)
            if ~isempty(chans2interpolate{c})
                % provide update
                fprintf('interpolating channel %s\n', chans2interpolate{c})
    
                % indentify the channel to interpolate
                chan_n = find(strcmp(params.labels, chans2interpolate{c}));
    
                % calculate distances with other electrodes
                chan_dist = -ones(length(dataset(1).processed(1).header.chanlocs), 1);
                for b = setdiff(1:chan_n, chan_n)
                    if dataset(1).processed(1).header.chanlocs(b).topo_enabled == 1
                        chan_dist(b) = sqrt((dataset(1).processed(1).header.chanlocs(b).X - dataset(1).processed(1).header.chanlocs(chan_n).X)^2 + ...
                            (dataset(1).processed(1).header.chanlocs(b).Y - dataset(1).processed(1).header.chanlocs(chan_n).Y)^2 + ...
                            (dataset(1).processed(1).header.chanlocs(b).Z - dataset(1).processed(1).header.chanlocs(chan_n).Z)^2);
                    end
                end
                chan_dist((chan_dist==-1)) = max(chan_dist);
                [~,chan_idx] = sort(chan_dist);
    
                % identify neighbouring channels
                chan_idx = chan_idx(1:params.interp_chans);
                chans2use = params.labels;
                chans2use = chans2use(chan_idx);
    
                % cycle through all datasets
                for d = 1:length(dataset(a).processed)
                    % select data
                    lwdata.header = dataset(a).processed(d).header;
                    lwdata.data = dataset(a).processed(d).data;
        
                    % interpolate using the neighboring electrodes
                    option = struct('channel_to_interpolate', chans2interpolate{c}, 'channels_for_interpolation_list', {chans2use}, ...
                        'suffix', '', 'is_save', 0);
                    lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
        
                    % update dataset
                    dataset(a).processed(d).header = lwdata.header;
                    dataset(a).processed(d).data = lwdata.data;  
                end
                
                % encode
                if c == 1
                    TEP_new(subject_idx).processing(9).process = sprintf('bad channels interpolated');
                    TEP_new(subject_idx).processing(9).date = sprintf('%s', date);
                end
                TEP_new(subject_idx).processing(9).params.bad{c} = chans2interpolate{c};
                TEP_new(subject_idx).processing(9).params.chans_used{c} = strjoin(chans2use, ' ');  
            end
        end
    else
        TEP_new(subject_idx).processing(9).process = sprintf('no channels interpolated');
        TEP_new(subject_idx).processing(9).date = sprintf('%s', date);
    end
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% re-reference and save as .set
fprintf('re-referencing and saving ...')
for a = 1:length(params.condition)
    for d = 1:length(dataset(a).processed)
        % select data
        lwdata.header = dataset(a).processed(d).header;
        lwdata.data = dataset(a).processed(d).data;

        % re-reference to common average
        option = struct('reference_list', {params.labels}, 'apply_list', {params.labels}, 'suffix', '', 'is_save', 0);
        lwdata = FLW_rereference.get_lwdata(lwdata, option);
        if a == 1 && d == 1
            TEP_new(subject_idx).processing(10).process = sprintf('re-referenced to common average');
            TEP_new(subject_idx).processing(10).date = sprintf('%s', date);
        end

        % export 
        export_EEGLAB(lwdata, lwdata.header.name, TEP_new(subject_idx).ID);
    end
end
fprintf('done.\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear a b c d prompt dlgtitle dims definput answer chans2interpolate chan_n chan_dist chan_idx chans2use lwdata option    
fprintf('section 3 finished.\n')

%% 4) SSP-SIR
% ----- section input -----
% -------------------------
fprintf('section 4: SSP-SIR\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear 
fprintf('section 4 finished.\n')

%% 5) ICA
% ----- section input -----
% -------------------------
fprintf('section 5: ICA\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear 
fprintf('section 5 finished.\n')

%% 6) group statistics
%% 7) visualization
%% functions
function dataset = reload_dataset(data2load, conditions, fieldname)
% =========================================================================
% Reloads pre-processed EEG data of a single subject for following 
% processing steps. 
% Input:    - list of datasets to load
%           - cell array with conditions
%           - fieldname
% =========================================================================  
% initiate output
dataset = struct;

% load all data
dpc = length(data2load)/2/length(conditions);
if mod(dpc, 1) == 0
    for c = 1:length(conditions)
        % note condition
        dataset(c).condition = conditions{c}; 

        % subset header and data files
        header_idx = logical([]);
        data_idx = logical([]);
        for d = 1:length(data2load)
            if contains(data2load(d).name, conditions{c}) 
                if contains(data2load(d).name, 'lw6') 
                    header_idx(d) = true;
                    data_idx(d) = false;
                elseif contains(data2load(d).name, 'mat') 
                    header_idx(d) = false;
                    data_idx(d) = true;
                end
            else
                header_idx(d) = false;
                data_idx(d) = false;
            end
        end
        headers = data2load(header_idx);
        datas = data2load(data_idx);

        % load all dataset for this condition
        if length(datas) == length(headers) && length(datas) == dpc
            for d = 1:dpc
                % load header
                load(sprintf('%s\\%s', headers(d).folder, headers(d).name), '-mat')
                statement = sprintf('dataset(c).%s(d).header = header;', fieldname);
                eval(statement) 
    
                % load data
                load(sprintf('%s\\%s', datas(d).folder, datas(d).name))
                statement = sprintf('dataset(c).%s(d).data = data;', fieldname);
                eval(statement) 
            end
        else
            error('ERROR: Wrong number of available datasets to load! Check manually.')
        end
    end
else
    error('ERROR: Wrong number of available datasets to load! Check manually.')
end
end
function export_EEGLAB(lwdata, filename, subj)
% =========================================================================
% exports data from letswave to EEGLAB format
% =========================================================================  
% dataset
EEG.setname = filename;
EEG.filename = [];
EEG.filepath = [];
EEG.subject = subj; 
EEG.session = 1;
    
% time properties
EEG.nbchan = lwdata.header.datasize(2);
EEG.trials = lwdata.header.datasize(1);
EEG.pnts = lwdata.header.datasize(6);
EEG.srate = 1/lwdata.header.xstep;
EEG.times = lwdata.header.xstart + (0:EEG.pnts-1)*lwdata.header.xstep;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end);
EEG.data = permute(single(lwdata.data),[2,6,1,3,4,5]);
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'SEEG_enabled');
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'topo_enabled');
    
% create events with appropriate latencies
EEG.event = lwdata.header.events;
if ~isempty(EEG.event)
    [EEG.event.type] = EEG.event.code;
    for e = 1:length(EEG.event)
        EEG.event(e).latency = (e-1)*EEG.pnts + 2001;
    end
    EEG.event = rmfield(EEG.event,'code');
end
    
% create required empty fields
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
save([filename,'.set'], 'EEG');
end