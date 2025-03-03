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
fprintf('done.\n\n')

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
        error('ERROR: This does not match with expected number of datasets (%d).\nPlease verify manually.\n', length(files2import))
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
fprintf('section 1 finished.\n\n')

%% 2) segment and pre-process for visual inspection
% ----- section input -----
params.prefix = 'dc ds crop';
params.suffix = {'no_mastoid' 'ep' 'dc' 'art_interp' 'preprocessed'};
params.eventcode = 'TMS';
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
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'raw');
    fprintf('done\n')
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
fprintf('done.\n\n')

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
        code = header.events(1).code;
        latency = header.events(1).latency;
        for c = 1:header.datasize(1)
            header.events(c).code = code;
            header.events(c).latency = latency;
            header.events(c).epoch = c;
        end
    
        % update dataset
        dataset(a).processed(b).header = header;
        dataset(a).processed(b).data = data;
    
        % save for letswave        
        save([header.name '.lw6'], 'header')
        save([header.name '.mat'], 'data')
    end
end
fprintf('done.\n\n')

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
fprintf('\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear data2load a b c d f lwdata option channel_all channel_mask channels2keep header data...
    code latency concat_idx fig_all open
fprintf('section 2 finished.\n\n')

%% 3) re-reference and export to EEGLAB
% ----- section input -----
params.prefix = 'preprocessed';
params.interp_chans = 6;
% -------------------------
fprintf('section 3: exporting to EEGLAB\n')

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('re-loading dataset... ')
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'processed');
    fprintf('done.\n')
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

% save as .set
fprintf('saving for EEGLAB... ')
for a = 1:length(params.condition)
    for d = 1:length(dataset(a).processed)
        fprintf('. ')

        % select data
        lwdata.header = dataset(a).processed(d).header;
        lwdata.data = dataset(a).processed(d).data;

        % export 
        export_EEGLAB(lwdata, lwdata.header.name, TEP_new(subject_idx).ID);
    end
end
fprintf('done.\n\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear a b c d e prompt dlgtitle dims definput answer chans2interpolate chan_n chan_dist chan_idx chans2use lwdata option 
fprintf('section 3 finished.\n\n')

%% 4) SSP-SIR
% ----- section input -----
params.prefix = 'preprocessed';
params.suffix = {'preprocessed' 'sspsir' 'ffilt'};
params.baseline = [-0.25 -0.006]; 
params.bandpass = [0.5, 80];
params.notch = [48.5, 51.5];
params.plot_output = true;
params.plot_toi = [-0.1 0.5];
% -------------------------
fprintf('section 4: SSP-SIR\n')

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('re-loading dataset... ')
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'processed');
    fprintf('done.\n')
end

% add eeglab to the top of search path and launch
addpath(fullfile(folder.toolbox, 'EEGLAB'));
eeglab

% apply SSP-SIR separately on data from each each session 
for a = 1:length(params.condition)
    fprintf('\n======================== %s session ========================\n', params.condition{a})

    % load all datasets, re-reference 
    fprintf('loading dataset: ')
    for b = 1:length(params.timepoint)
        fprintf('. ')
        % load dataset, re-reference  
        name = sprintf('%s %s %s %s %s.set', params.suffix{1}, study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b}); 
        EEG = pop_loadset('filename', name, 'filepath', folder.processed);
        EEG.filename = name;
        EEG.filepath = folder.processed;
        EEG = pop_reref(EEG, []);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, b);
        eeglab redraw  
        if a == 1 && b == 1
            TEP_new(subject_idx).processing(10).process = 're-referenced to common average';
            TEP_new(subject_idx).processing(10).date = sprintf('%s', date);
        end
    end
    fprintf('done\n')

    % merge the data
    fprintf('merging ...\n')
    merged_EEG = pop_mergeset(ALLEEG, 1:length(ALLEEG), 0);

    % apply SSP-SIR - spherical model 
    fprintf('applying SSP-SIR ...\n')
    merged_EEG = pop_tesa_sspsir(merged_EEG, 'artScale', 'automatic', 'PC', []);
    prompt = {'number of rejected PCs:'};
    dlgtitle = 'SSP-SIR';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    if a == 1 
        TEP_new(subject_idx).processing(11).process = 'muscular artifact removed using SSP-SIR';
        TEP_new(subject_idx).processing(11).params.method = 'SSP-SIR';
        TEP_new(subject_idx).processing(11).params.leadfield = 'default spherical';
        TEP_new(subject_idx).processing(11).suffix = params.suffix{2};
        TEP_new(subject_idx).processing(11).date = sprintf('%s', date);
    end
    TEP_new(subject_idx).processing(11).params.PC(a).condition = params.condition{a};
    TEP_new(subject_idx).processing(11).params.PC(a).removed = str2num(input{1,1});

    % split back to original datasets
    idx_start = 1;
    for b = 1:length(params.timepoint)
        % identify number of epochs 
        n_epochs = length(ALLEEG(b).event);

        % extract epochs and update ALLEEG
        EEG = pop_select(merged_EEG, 'trial', idx_start:(idx_start + n_epochs - 1));

        % save new dataset
        name = sprintf('%s %s %s %s %s', params.suffix{2}, study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b}); 
        EEG.setname = name;
        EEG.filename = sprintf('%s.set', name);
        pop_saveset(EEG, 'filename', name, 'filepath', folder.processed);

        % update starting index 
        idx_start = idx_start + n_epochs;       
    end
    
    % apply frequency filters and save
    fprintf('final pre-processing:\n')
    for b = 1:length(params.timepoint)
        fprintf('%s - %s:\n', params.condition{a}, params.timepoint{b})

        % load dataset
        name = sprintf('%s %s %s %s %s.set', params.suffix{2}, study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b}); 
        EEG = pop_loadset('filename', name, 'filepath', folder.processed);
        [ALLEEG, EEG, b] = eeg_store(ALLEEG, EEG, b);
        eeglab redraw 

        % bandpass filter
        EEG = pop_eegfiltnew(EEG, 'locutoff', params.bandpass(1), 'hicutoff', params.bandpass(2));

        % notch filter
        EEG = pop_eegfiltnew(EEG, 'locutoff', params.notch(1), 'hicutoff', params.notch(2), 'revfilt', 1);

        % % check spectrum
        % pop_spectopo(EEG, 1, [], 'EEG', 'freqrange', [0 100]);

        % encode to info structure
        if a == 1 && b == 1
            TEP_new(subject_idx).processing(12).process = sprintf('frequency filters applied');
            TEP_new(subject_idx).processing(12).params.bandpass = params.bandpass;
            TEP_new(subject_idx).processing(12).params.notch = params.notch;
            TEP_new(subject_idx).processing(12).suffix = params.suffix{3};
            TEP_new(subject_idx).processing(12).date = sprintf('%s', date);
        end

        % update dataset and save for letswave
        name = sprintf('%s %s %s %s %s %s %s %s', params.suffix{3}, params.suffix{2},...
            study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b}); 
        lwdata = export_lw(EEG, dataset(a).processed(b).header, name);
        dataset(a).sspsir(b).header = lwdata.header;
        dataset(a).sspsir(b).data = lwdata.data;
        header = lwdata.header;
        data = lwdata.data;
        save([name '.lw6'], 'header')
        save([name '.mat'], 'data')
        fprintf('\n')
    end

    % clear data strucure for next iteration
    ALLEEG = []; 
end

% update figure counter
fig_all = findall(0, 'Type', 'figure');
figure_counter = length(fig_all) + 1;

% plot the output figure
if params.plot_output
    % define common visual parameters
    params.labels = {dataset(1).sspsir(1).header.chanlocs.labels};
    visual.x = dataset(1).sspsir(1).header.xstart : dataset(1).sspsir(1).header.xstep : ...
        dataset(1).sspsir(1).header.xstep * dataset(1).sspsir(1).header.datasize(6) + dataset(1).sspsir(1).header.xstart - dataset(1).sspsir(1).header.xstep;
    visual.labels = {'original', 'filtered'};
    visual.xlim = params.plot_toi;
    visual.colors = [0.3333    0.4471    0.9020;
                    1.0000    0.0745    0.6510];

    % launch the figure
    fig = figure(figure_counter);
    set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1])
    hold on

    % plot data
    for a = 1:length(params.condition)
        % identify EOI
        if strcmp(TEP_new(subject_idx).stimulation(a).hemisphere, 'right')
            visual.text  = 'C4';
        elseif strcmp(TEP_new(subject_idx).stimulation(a).hemisphere, 'left')
            visual.text  = 'C3';
        end
        eoi = find(strcmp(params.labels, visual.text));

        % cycle through timepoints
        for b = 1:length(params.timepoint)
            % define t value
            visual.t_value = tinv(0.975, size(dataset(a).processed(b).data, 1) - 1); 

            % select original data
            visual.data(1, :) = squeeze(mean(dataset(a).processed(b).data(:,eoi,1,1,1,:), 1));  
            visual.sem(1, :) = squeeze(std(dataset(a).processed(b).data(:,eoi,1,1,1,:), 0, 1)) / sqrt(size(dataset(a).processed(b).data, 1)); 
            visual.CI_upper(1, :) = visual.data(1, :) + visual.t_value * visual.sem(1, :); 
            visual.CI_lower(1, :) = visual.data(1, :) - visual.t_value * visual.sem(1, :); 
            
            % select filtered data
            visual.data(2, :) = squeeze(mean(dataset(a).sspsir(b).data(:,eoi,1,1,1,:), 1));  
            visual.sem(2, :) = squeeze(std(dataset(a).sspsir(b).data(:,eoi,1,1,1,:), 0, 1)) / sqrt(size(dataset(a).sspsir(b).data, 1)); 
            visual.CI_upper(2, :) = visual.data(2, :) + visual.t_value * visual.sem(2, :); 
            visual.CI_lower(2, :) = visual.data(2, :) - visual.t_value * visual.sem(2, :); 
            
            % plot
            subplot(length(params.condition), length(params.timepoint), (a-1)*length(params.timepoint) + b)
            plot_ERP(visual, 'xlim', visual.xlim, 'colours', visual.colors, 'labels', visual.labels)
            title(sprintf('%s - %s', params.condition{a}, params.timepoint{b}))
        end
    end

    % save and update figure counter
    saveas(fig, sprintf('%s\\figures\\SSPSIR_%s.png', folder.output, TEP_new(subject_idx).ID))
    figure_counter = figure_counter + 1;
end

% open letswave if not already open
addpath(genpath([folder.toolbox '\letswave 6']));
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
fprintf('\n')

% save and continue
save(output_file, 'TEP_new','-append')
clear a b c f e fig_all name match prompt discarded answer data header lwdata data2load definput dims dlgtitle eoi fig idx_start input latency code ...
    n_epochs open tmpEEG tmpstr visual ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG merged_EEG globalvars LASTCOM PLUGINLIST STUDY GUI_handle GUI_OK
fprintf('section 4 finished.\nPlease check for bad trials now.\n\n')

%% 5) ICA
% ----- section input -----
params.prefix = 'ar ffilt sspsir';
params.suffix = {'ar' 'ica'};
params.ICA_comp = 25;
% -------------------------
fprintf('section 5: ICA\n')

% load dataset with bad trials removed
fprintf('updating dataset... ')
if exist('dataset') ~= 1
    % load previous dataset
    data2load = dir(sprintf('%s*%s*', params.prefix(4:end), TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'sspsir');

    % append new dataset
    dataset_old = dataset;
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'checked');
    for a = 1:length(params.condition)
        dataset(a).sspsir = dataset_old(a).sspsir;
    end
    clear dataset_old
else
    % append new dataset
    dataset_old = dataset;
    data2load = dir(sprintf('%s*%s*', params.prefix, TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'checked');
    for a = 1:length(params.condition)
        dataset(a).sspsir = dataset_old(a).sspsir;
    end
    clear dataset_old
end
fprintf('done.\n')

% encode bad trials
fprintf('encoding bad trials... ')
TEP_new(subject_idx).processing(13).process = sprintf('bad trials discarded');
TEP_new(subject_idx).processing(13).params.GUI = 'letswave';
TEP_new(subject_idx).processing(13).suffix = params.suffix{1};
TEP_new(subject_idx).processing(13).date = sprintf('%s', date);
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        % load header
        load(sprintf('%s %s %s %s %s .lw6', params.prefix, study, TEP_new(subject_idx).ID, params.condition{a}, params.timepoint{b}), '-mat') 

        % extract discarded expochs
        if ~isempty(header.history(end).configuration)
            if ~isempty(header.history(end).configuration.parameters.rejected_epochs)
                discarded = header.history(end).configuration.parameters.rejected_epochs;
            end
        end

        % encode 
        TEP_new(subject_idx).processing(13).params.discarded{a, b} = discarded;
        TEP_new(subject_idx).processing(13).params.kept(a, b) = size(dataset(a).checked(b).data, 1);
    end
end
fprintf('done\n')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% compute ICA matrix and save for letswave
for a = 1:length(params.condition)
    fprintf('\n======================== %s session ========================\n', params.condition{a})

    % select dataset
    lwdataset = dataset(a).checked;

    % compute ICA and save  
    fprintf('computing ICA matrix:\n')
    option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix{2}, 'is_save', 1);
    lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    fprintf('done\n')

    % extract ICA parameters
    matrix(a).condition = params.condition{a};
    matrix(a).mix = lwdataset(1).header.history(end).option.mix_matrix;
    matrix(a).unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
    if a == 1
        params.ICA_chanlocs = lwdataset(1).header.chanlocs;
        for i = 1:size(matrix(a).mix, 2)
            params.ICA_labels{i} = ['IC',num2str(i)];
        end
        params.ICA_SR = 1/lwdataset(1).header.xstep;
    end

    % update dataset and adjust for letswave 6
    dataset(a).ica = lwdataset;
    for b = 1:length(params.timepoint)
        dataset(a).ica(b).header.history(10).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
        dataset(a).ica(b).header.history(10).configuration.parameters = dataset(a).ica(b).header.history(10).option;  
        [dataset(a).ica(b).header.history(10).configuration.parameters.ICA_um] = dataset(a).ica(b).header.history(10).configuration.parameters.unmix_matrix; 
        [dataset(a).ica(b).header.history(10).configuration.parameters.ICA_mm] = dataset(a).ica(b).header.history(10).configuration.parameters.mix_matrix; 
        dataset(a).ica(b).header.history(10).configuration.parameters = rmfield(dataset(a).ica(b).header.history(10).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
        header = dataset(a).ica(b).header;
        save(sprintf('%s.lw6', dataset(a).ica(b).header.name), 'header');
    end

    % unmix data
    for b = 1:length(dataset(a).ica)
        for e = 1:size(dataset(a).ica(b).data, 1)
            dataset(a).unmixed(b).header = dataset(a).ica(b).header;
            dataset(a).unmixed(b).data(e, :, 1, 1, 1, :) = matrix(a).unmix * squeeze(dataset(a).ica(b).data(e, :, 1, 1, 1, :));        
        end
    end
end

% update info structure
TEP_new(subject_idx).processing(14).process = 'ICA matrix computed';
TEP_new(subject_idx).processing(14).params.method = 'runica';
TEP_new(subject_idx).processing(14).params.components = params.ICA_comp;
TEP_new(subject_idx).processing(14).params.chanlocs = params.ICA_chanlocs;
TEP_new(subject_idx).processing(14).params.labels = params.ICA_labels;
TEP_new(subject_idx).processing(14).params.SR = params.ICA_SR;
TEP_new(subject_idx).processing(14).params.matrix = matrix;
TEP_new(subject_idx).processing(14).suffix = params.suffix{2};
TEP_new(subject_idx).processing(14).date = sprintf('%s', date);

% plot IC spectral content separately for each session
fprintf('estimating spectral content...\n')
for a = 1:length(params.condition)
    % calculate PSD across all timepoints, components and trials 
    psd = [];
    for b = 1:length(params.timepoint)
        for c = 1:params.ICA_comp
            for e = 1:size(dataset(a).unmixed(b).data, 1)
                [psd(b, c, e, :), freq] = pwelch(squeeze(dataset(a).unmixed(b).data(e, c, 1, 1, 1, :)), ...
                    [], [], [], TEP_new(subject_idx).processing(14).params.SR);  
            end
        end
    end
    TEP_new(subject_idx).processing(14).params.PSD = squeeze(mean(psd, [1, 3]));

    % plot component topographies and spectral content
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    for f = 1:params.ICA_comp
        % plot the topography
        matrix = TEP_new(subject_idx).processing(14).params.matrix(a).mix;
        labels = {dataset(1).ica(1).header.chanlocs.labels};
        subplot(ceil(params.ICA_comp/3), 6, (f-1)*2 + 1);
        topoplot(double(matrix(:, f)'), params.ICA_chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(params.ICA_labels{f})
    
        % plot the psd
        subplot(ceil(params.ICA_comp/3), 6, (f-1)*2 + 2);
        plot(freq(1:82), log10(TEP_new(subject_idx).processing(14).params.PSD(f, 1:82)));
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
    end
    sgtitle(sprintf('%s - %s session', TEP_new(subject_idx).ID, params.condition{a}))
    saveas(gcf, sprintf('%s\\figures\\ICA_%s_%s.png', folder.output, TEP_new(subject_idx).ID, params.condition{a}));
end
fprintf('done.\n')

% open letswave 6 
addpath(genpath([folder.toolbox '\letswave 6']));
letswave

% save and continue
save(output_file, 'TEP_new','-append')
clear a b c d e f i data2load check discarded match prompt definput input input_old dims dlgtitle matrix fig_all open ...
    psd freq option lwdataset labels header
fprintf('section 5 finished.\nPlease procede to ICA now\n\n')

%% 6) encode ICA
% ----- section input -----
params.suffix = {'icfilt'};
params.ICA_comp = 25;
params.plot_toi = [-0.1 0.5];
params.eoi = 'Cz';
% -------------------------
fprintf('section 6: encode ICA\n')

% load dataset if necessary
if exist('dataset') ~= 1
    fprintf('loading dataset... ')
    data2load = dir(sprintf('ica*%s*', TEP_new(subject_idx).ID));
    dataset = reload_dataset(data2load, params.condition, 'ica');
    fprintf('done\n')
end

% encode 
TEP_new(subject_idx).processing(15).process = 'artifactual ICs discarded';
TEP_new(subject_idx).processing(15).suffix = params.suffix{1};
TEP_new(subject_idx).processing(15).date = sprintf('%s', date);
for a = 1:length(params.condition)
    % ask for the input
    prompt = {'blinks:', 'horizontal:', 'TMS:', 'muscles:', 'slow artifacts:'};
    dlgtitle = sprintf('ICA - %s session', params.condition{a});  
    dims = [1 60];
    definput = {'', '', '', '', ''};
    input = inputdlg(prompt,dlgtitle,dims,definput);

    % encode
    TEP_new(subject_idx).processing(15).params.kept(a).condition = params.condition{a};
    TEP_new(subject_idx).processing(15).params.kept(a).components = params.ICA_comp - length([str2num(input{1}), str2num(input{2}), str2num(input{3}), str2num(input{4}), str2num(input{5})]);
    TEP_new(subject_idx).processing(15).params.removed(a).condition = params.condition{a};
    TEP_new(subject_idx).processing(15).params.removed(a).blinks = str2num(input{1});
    TEP_new(subject_idx).processing(15).params.removed(a).horizontal = str2num(input{2});
    TEP_new(subject_idx).processing(15).params.removed(a).TMS = str2num(input{3});
    TEP_new(subject_idx).processing(15).params.removed(a).muscles = str2num(input{4});
    TEP_new(subject_idx).processing(15).params.removed(a).slow = str2num(input{5});
end
save(output_file, 'TEP_new', '-append')

% launch the figure
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3) / 2.5, screen_size(4) / 2])

% define common visual parameters
params.labels = {dataset(1).ica(1).header.chanlocs.labels};
visual.x = dataset(1).ica(1).header.xstart : dataset(1).ica(1).header.xstep : ...
    dataset(1).ica(1).header.xstep * dataset(1).ica(1).header.datasize(6) + dataset(1).ica(1).header.xstart - dataset(1).ica(1).header.xstep;
visual.labels = params.condition;
visual.xlim = params.plot_toi;
visual.colors = [0.9216    0.1490    0.1490;
                0    0.4471    0.7412];

% extract mean signal per condition
for a = 1:length(params.condition)
    % identify EOI
    if strcmp(params.eoi, 'target')
        if strcmp(TEP_new(subject_idx).stimulation(a).hemisphere, 'right')
            visual.text  = 'C4';
        elseif strcmp(TEP_new(subject_idx).stimulation(a).hemisphere, 'left')
            visual.text  = 'C3';
        end        
    else
        visual.text  = params.eoi;
    end
    eoi = find(strcmp(params.labels, visual.text));

    % subset data
    data2plot = [];
    for b = 1:length(params.timepoint)
        load(sprintf('%s %s.mat', params.suffix{1}, dataset(a).ica(b).header.name))
        data2plot(end + 1 : end + size(data,1), :) = squeeze(data(:, eoi, 1, 1, 1, :));
    end

    % define t value
    visual.t_value = tinv(0.975, size(data2plot, 1) - 1); 

    % select original data
    visual.data(a, :) = mean(data2plot, 1)';  
    visual.sem(a, :) = std(data2plot, 0, 1)' / sqrt(size(data2plot, 1)); 
    visual.CI_upper(a, :) = visual.data(a, :) + visual.t_value * visual.sem(a, :); 
    visual.CI_lower(a, :) = visual.data(a, :) - visual.t_value * visual.sem(a, :);
end

% plot
plot_ERP(visual, 'xlim', visual.xlim, 'colours', visual.colors, 'labels', visual.labels)

% save and update counter
saveas(fig, sprintf('%s\\figures\\TEP_avg_%s.png', folder.output, TEP_new(subject_idx).ID))
figure_counter = figure_counter + 1;

% ask for continuation
fprintf('section 6 finished.\n')
answer = questdlg('Do you want to continue with next subject?', 'Continue?', 'YES', 'NO', 'YES'); 
if strcmp(answer, 'YES')
    subject_idx = subject_idx + 1;
    clear dataset
    fprintf('proceeding to subject %d.\n\n', subject_idx)
end
clear a b data2load prompt definput input dims dlgtitle answer fig screen_size visual eoi data data2plot

%% 7) load and plot group data
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.subjects = 20;
params.baseline = [-0.3 -0.005];
% -------------------------
fprintf('section 7: load and plot group data\n')

% update
% load(output_file, 'TEP_new_data')
% dataset = TEP_new_data;
clear dataset subject_idx
cd(folder.output)

% load data of all subjects
fprintf('loading data: ')
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        fprintf('. ')
        for s = 1:params.subjects            
            % load the data
            load(sprintf('%s %s %s %s %s .mat', params.prefix, study, TEP_new(s).ID, params.condition{a}, params.timepoint{b}))

            % append to the dataset
            dataset.avg(a, b, s, :, :) = squeeze(mean(data, 1));
        end
    end
end
fprintf('done.\n')

% load dataset information
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')
params.labels = {header.chanlocs.labels};
params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
params.chanlocs = header.chanlocs;

% perform baseline normalization to z-scores
fprintf('correcting for baseline: ')
baseline.idx = find(params.x >= params.baseline(1) & params.x <= params.baseline(2));
dataset.normalized = dataset.avg;
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        fprintf('. ')
        for s = 1:params.subjects 
            for c = 1:length(params.chanlocs)
                % extract baseline segment
                baseline.values = squeeze(dataset.normalized(a, b, s, c, baseline.idx));

                % compute mean and standard deviation of baseline
                baseline.mean = mean(baseline.values, 'all');
                baseline.std = std(baseline.values, 0, 'all');

                % normalize the signal
                if baseline.std > 0  
                    dataset.normalized(a, b, s, c, :) = (squeeze(dataset.normalized(a, b, s, c, :)) - baseline.mean) / baseline.std;
                else
                    warning('Baseline std = 0 for %s condition, timepoint %s, subject %d, electrode %s. No normalization applied.', ...
                        params.condition{a}, params.timepoint{b}, s, params.labels{c});
                end
                
            end
        end
    end
end
fprintf('done.\n')

% prepare flip dictionary
labels_flipped = params.labels;
for i = 1:length(params.labels)
    electrode_n = str2num(params.labels{i}(end));
    if isempty(electrode_n)
    else
        if mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
            a = find(strcmpi(params.labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        else                                    % even number --> right hemisphere 
            label_new = params.labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
            a = find(strcmpi(params.labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        end
    end
end
labels_dict = cat(1, params.labels, labels_flipped)';

% flip normalized data to homogenize side of stimulation --> if right, flip
fprintf('flipping data: ')
addpath(genpath([folder.toolbox '\letswave 6']));
header.datasize(1) = 1;
for s = 1:params.subjects 
    fprintf('. ')
    for a = 1:length(params.condition)
        for b = 1:length(params.timepoint)
            if contains(TEP_new(s).stimulation(a).hemisphere, 'right')
                data(1, :, 1, 1, 1, :) = squeeze(dataset.normalized(a, b, s, :, :));
                [header, data, ~] = RLW_flip_electrodes(header, data, labels_dict);
                dataset.normalized_flipped(a, b, s, :, :) = squeeze(data);
            else
                dataset.normalized_flipped(a, b, s, :, :) = dataset.normalized(a, b, s, :, :);
            end
        end
    end
end
fprintf('done.\n')

% % export to letswave to check the flip
% header.name = 'original';
% data(1, :, 1, 1, 1, :) = squeeze(dataset.normalized(2, 1, 1, :, :));
% save('original.mat', 'data'); save('original.lw6', 'header'); 
% header.name = 'flipped';
% data(1, :, 1, 1, 1, :) = squeeze(dataset.normalized_flipped(2, 1, 1, :, :));
% save('flipped.mat', 'data'); save('flipped.lw6', 'header'); 
% letswave

% calculate change from baseline + its GFP
for a = 1:length(params.condition)
    for b = 2:length(params.timepoint)
        for s = 1:params.subjects
            dataset.change(a, b-1, s, :, :) = dataset.normalized_flipped(a, b, s, :, :) - dataset.normalized_flipped(a, 1, s, :, :);
            dataset.change_GFP(a, b-1, s, :) = std(dataset.change(a, b-1, s, :, :), 1, 4);
        end
    end
end

% define common visual parameters
visual.x = params.x;
visual.labels = params.condition;
visual.chanlabels = params.labels;
visual.t_value = tinv(0.975, size(dataset.change_GFP, 3) - 1); 

% ============== overall average TEP ==============
fprintf('plotting overall average TEP... \n')

% extract overall mean values
for c = 1:size(dataset.normalized_flipped, 4)
    dataset.visual_mean.data(c, :) = squeeze(mean(dataset.normalized_flipped(:, :, :, c, :), 1:3))'; 
    dataset.visual_mean.sem(c, :) = squeeze(std(dataset.normalized_flipped(:, :, :, c, :), 0, 1:3))'; 
    dataset.visual_mean.CI_upper(c, :) = dataset.visual_mean.data(c, :) + visual.t_value * dataset.visual_mean.sem(c, :); 
    dataset.visual_mean.CI_lower(c, :) = dataset.visual_mean.data(c, :) - visual.t_value * dataset.visual_mean.sem(c, :);
end

% launch the figure
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 2])

% select plotting colours
for c = 1:length(params.chanlocs)
    visual.colors(c, :) = [0.3020    0.7451    0.9333];
end

% select the data
visual.data = dataset.visual_mean.data;
visual.CI_upper = dataset.visual_mean.CI_upper;
visual.CI_lower = dataset.visual_mean.CI_lower;
    
% plot 
plot_ERP(visual, 'xlim', [-0.1 0.4], 'colours', visual.colors, 'legend', 'off', 'shading', 'off', 'interpolated', [-0.005 0.01], 'eoi', 'Cz')

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\group_avg.svg', folder.output))
figure_counter = figure_counter + 1;

% ============== overall average GFP ==============
fprintf('plotting overall average GFP... \n')

% launch the figure
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 2])

% select the data
visual.data = std(dataset.visual_mean.data, 1, 1); 
    
% plot 
plot_ERP(visual, 'xlim', [-0.1 0.4], 'colours', [0 0 0], 'legend', 'off', 'shading', 'off', 'interpolated', [-0.005 0.01])

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\group_GFP.svg', folder.output))
figure_counter = figure_counter + 1;

% ============ plot change - butterfly ============
fprintf('plotting change... \n')

% extract mean change per condition, timepoint, and electrode
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)-1   
        for c = 1:length(params.chanlocs)
            dataset.visual_butterfly.data(a, b, c, :) = squeeze(mean(dataset.change(a, b, :, c, :), 3))';  
            dataset.visual_butterfly.sem(a, b, c, :) = squeeze(std(dataset.change(a, b, :, c, :), 0, 3))' / sqrt(size(dataset.change(a, b, :, c, :), 3)); 
            dataset.visual_butterfly.CI_upper(a, b, c, :) = dataset.visual_butterfly.data(a, b, c, :) + visual.t_value * dataset.visual_butterfly.sem(a, b, c, :); 
            dataset.visual_butterfly.CI_lower(a, b, c, :) = dataset.visual_butterfly.data(a, b, c, :) - visual.t_value * dataset.visual_butterfly.sem(a, b, c, :);
        end
    end 
end

% launch the figure
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [0, screen_size(4)/8, screen_size(3), screen_size(4) / 1.5])

% plot mean change from baseline
visual.ylim = [];
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)-1 
        % select plotting colours
        for c = 1:length(params.chanlocs)
            if a == 1
                visual.colors(c, :) = [0.9216    0.1490    0.1490];
            elseif a == 2
                visual.colors(c, :) = [0    0.4471    0.7412];
            end
        end

        % select the data
        visual.data = squeeze(dataset.visual_butterfly.data(a, b, :, :));
        visual.CI_upper = squeeze(dataset.visual_butterfly.CI_upper(a, b, :, :));
        visual.CI_lower = squeeze(dataset.visual_butterfly.CI_lower(a, b, :, :));
    
        % plot averages of both conditions in one plot
        subplot(2, length(params.timepoint)-1, (a-1)*(length(params.timepoint)-1) + b)
        plot_ERP(visual, 'xlim', [-0.025 0.1], 'colours', visual.colors, 'legend', 'off', 'shading', 'off')
        if a == 1
            title(sprintf('%s', params.timepoint{b+1}))
        end
        hold on

        % get y limits
        visual.ylim(a, b, :) = get(gca, 'ylim');
    end
end

% regulate the resolution
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)-1
        subplot(2, length(params.timepoint)-1, (a-1)*(length(params.timepoint)-1) + b)
        ylim([min(visual.ylim(:, :, 1), [], 'all'), max(visual.ylim(:, :, 2), [], 'all')])
        line([0, 0], [min(visual.ylim(:, :, 1), [], 'all'), max(visual.ylim(:, :, 2), [], 'all')], 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')
        legend('off');
    end
end

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\group_change_butterfly.png', folder.output))
figure_counter = figure_counter + 1;

% ================ plot change GFP ================
fprintf('plotting GFP of change... \n')

% extract mean GFP change per condition and timepoint
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)-1         
        dataset.visual_GFP.data(a, b, :) = squeeze(mean(dataset.change_GFP(a, b, :, :), 3))';  
        dataset.visual_GFP.sem(a, b, :) = squeeze(std(dataset.change_GFP(a, b, :, :), 0, 3))' / sqrt(size(dataset.change_GFP(a, b, :, :), 3)); 
        dataset.visual_GFP.CI_upper(a, b, :) = dataset.visual_GFP.data(a, b, :) + visual.t_value * dataset.visual_GFP.sem(a, b, :); 
        dataset.visual_GFP.CI_lower(a, b, :) = dataset.visual_GFP.data(a, b, :) - visual.t_value * dataset.visual_GFP.sem(a, b, :);
    end 
end

% launch the figure
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [0, screen_size(4)/4, screen_size(3), screen_size(4) / 2])

% select plotting colours
visual.colors = [0.9216    0.1490    0.1490;
                0    0.4471    0.7412];

% plot mean GFP change from baseline
visual.ylim = [];
for b = 1:length(params.timepoint)-1 
    % select the data
    visual.data = squeeze(dataset.visual_GFP.data(:, b, :));
    visual.CI_upper = squeeze(dataset.visual_GFP.CI_upper(:, b, :));
    visual.CI_lower = squeeze(dataset.visual_GFP.CI_lower(:, b, :));

    % plot averages of both conditions in one plot
    subplot(1, length(params.timepoint)-1, b)
    plot_ERP(visual, 'xlim', [-0.025 0.1], 'colours', visual.colors, 'labels', visual.labels)
    title(sprintf('%s', params.timepoint{b+1}))
    hold on

    % get y limits
    visual.ylim(b, :) = get(gca, 'ylim');
end

% regulate the resolution
for b = 1:length(params.timepoint)-1
    subplot(1, length(params.timepoint)-1, b)
    ylim([0, max(visual.ylim(:, 2))])
    line([0, 0], [0, max(visual.ylim(:, 2))], 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')
    legend(findobj(gca, 'Type', 'line', '-not', 'Color', 'black'));
end

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\group_change_GFP.png', folder.output))
figure_counter = figure_counter + 1;

% save the dataset and clean up
TEP_new_data = dataset;
save(output_file, 'TEP_new_data', '-append')
clear a b c s i data header data baseline fig screen_size visual electrode_n ...
    labels_flipped labels_dict label_new
fprintf('section 7 finished.\n')

%% 8) export for Ragu
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.subjects = 20;
params.toi = [-0.1 0.4];
% -------------------------
fprintf('section 8: export for Ragu\n')

% update
load(output_file, 'TEP_new_data')
dataset = TEP_new_data;

% calculate crop limits
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')                   
x_start = (params.toi(1) - header.xstart)/header.xstep + 1;
x_end = (params.toi(2) - header.xstart)/header.xstep;

% write text files for Ragu
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)  
        for s = 1:params.subjects
            % select data
            data = squeeze(dataset.normalized_flipped(a, b, s, :, x_start:x_end))';

            % save as .csv               
            name = sprintf('%s_%s_%s_%s.csv', study, TEP_new(s).ID, params.condition{a}, params.timepoint{b}); 
            writematrix(data, sprintf('%s\\export\\%s', folder.output, name))
        end
    end
end

% create the montage file
name = sprintf('%s\\export\\%s_montage.xyz', folder.output, study);  
fileID = fopen(name, 'a');
fprintf(fileID, '30\r\n');
for c = 1:length(params.chanlocs)
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        params.chanlocs(c).X, params.chanlocs(c).Y, params.chanlocs(c).Z, params.chanlocs(c).labels);
end
fclose(fileID)

% clear and move on
clear a b c s header x_start x_end data name fileID
fprintf('section 8 finished.\n')

%% 9) import microstates
% ----- section input -----
params.filenames = {'CAPSTEP_new.mat' 'microstates_export.xlsx'};
params.MS_labels = {'N15' 'P180' 'P30' 'N45_P60' 'N100' '' ''};
params.MS_variables = {'On' 'Off' 'Dur' 'AUC' 'COG' 'GFP'};
params.MS_measures = {'onset' 'offset' 'duration' 'AUC' 'COG' 'GFP'};
params.MS_manual.N15 = [11, 21];
params.MS_manual.P180 = [143, 234];
params.MS_manual.P30 = [25, 40];
params.MS_manual.N45_P60 = [41, 84];
params.MS_manual.N100 = [85, 122];
% -------------------------
fprintf('section 9: import microstate info\n')

% load & save microstate maps
load(sprintf('%s\\microstates\\%s', folder.output, params.filenames{1}))
for a = 1:length(params.MS_labels)
    if isempty(params.MS_labels{a})
        idx(a) = false;
    else
        idx(a) = true;
    end
end
TEP_new_measures.microstates.labels = params.MS_labels(idx);
TEP_new_measures.microstates.maps = rd.MSMaps(idx, :);
save(output_file, 'TEP_new_measures', '-append')

% identify export Excel
export.file = sprintf('%s\\microstates\\%s', folder.output, params.filenames{2});
[~, export.sheets] = xlsfinfo(export.file);

% identify conditions
counter = 1;
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint) 
        params.MS_conditions{counter} = sprintf('%s_%s', params.condition{a}, params.timepoint{b});
        counter = counter + 1;
    end
end
TEP_new_measures.microstates.conditions = params.MS_conditions;

% initiate output fields
for c = 1:length(params.MS_labels) 
    if idx(c)
        for d = 1:length(params.MS_measures)
            statement = sprintf('TEP_new_measures.microstates.%s.%s = [];', params.MS_labels{c}, params.MS_measures{d});
            eval(statement)
        end
    end
end

% extract data from the Excel
for e = 1:numel(export.sheets)
    if contains(export.sheets{e}, 'Class')
        % identify microstate class
        export.class = str2double(regexp(export.sheets{e}, '\d+', 'match')); 
        if idx(export.class)
            % identify microstate label
            export.label = params.MS_labels{export.class};
    
            % read the sheet as a table
            export.table = readtable(export.file, 'Sheet', export.sheets{e});

            % extract all variables
            for c = 1:length(params.MS_measures)
                for d = 1:length(params.MS_conditions)
                    % extract the right column
                    export.column = ismember(export.table.Properties.VariableNames, sprintf('%s_%s', params.MS_conditions{d}, params.MS_variables{c}));
                    export.values = table2cell(export.table(:, export.column));
                    
                    % save to the output structure
                    for s = 1:height(export.table)                        
                        statement = sprintf('TEP_new_measures.microstates.%s.%s(s).%s = export.values{s};', ...
                            params.MS_labels{export.class}, params.MS_measures{c}, params.MS_conditions{d});
                        eval(statement)
                    end
                end
            end

            % calculate average across all subjects and conditions
            for c = 1:length(params.MS_measures)
                % calculate the mean
                data_all = [];
                for d = 1:length(params.MS_conditions)
                    statement = sprintf('data = [TEP_new_measures.microstates.%s.%s.%s];', ...
                        params.MS_labels{export.class}, params.MS_measures{c}, params.MS_conditions{d});
                    eval(statement)
                    data_all = [data_all data];
                end
                
                % save to the output structure
                statement = sprintf('TEP_new_measures.microstates.%s.mean.%s = nanmean(data);', ...
                    params.MS_labels{export.class}, params.MS_measures{c}, params.MS_labels{export.class}, params.MS_measures{c});
                eval(statement)
            end
        end
    end
end

% encode manually measured time windows
for a = 1:length(params.MS_labels)
    if idx(a)
        statement = sprintf('TEP_new_measures.microstates.%s.toi = params.MS_manual.%s;', params.MS_labels{a}, params.MS_labels{a});
        eval(statement)
    end
end
save(output_file, 'TEP_new_measures', '-append')

% clear and move on
clear a b c d e s rd idx export counter statement data data_all
fprintf('section 9 finished.\n')

%% 10a) statistics: change from baseline - separately per electrode
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.toi = [0.01 0.4];
params.permutations = 1000;
params.alpha = 0.05;
params.plot = [-0.025 0.4];
params.components = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
params.N45_offset = 58;
params.topo_lim = [-2, 2];
% -------------------------
fprintf('section 10: statistics - change from baseline\n')

% update
load(output_file, 'TEP_new_data', 'TEP_new_measures')
dataset = TEP_new_data.normalized_flipped;
cd(folder.output)
clear dataset subject_idx

% load dataset information
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')
params.labels = {header.chanlocs.labels};
params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
params.chanlocs = header.chanlocs;

% identify analysed time window
toi = find(params.x >= params.toi(1) & params.x <= params.toi(2));

% initialize output
p_values = NaN(size(dataset, 1), size(dataset, 2) - 1, size(dataset, 4), length(toi));
p_values_correct = ones(size(p_values));
t_values = NaN(size(p_values));

% perform sample-wise paired t-test
fprintf('computing sample-wise paired t-test: ')
for a = 1:size(dataset, 1)
    for b = 2:size(dataset, 2)
        fprintf('. ')
        for c = 1:size(dataset, 4)
            for d = 1:length(toi)
                % extract pre-patch and post-patch values across subjects
                pre_patch = squeeze(dataset(a, 1, :, c, toi(d))); 
                post_patch = squeeze(dataset(a, b, :, c, toi(d))); 
                
                % perform paired t-test
                [~, p, ~, stats] = ttest(post_patch, pre_patch);
                
                % store results
                p_values(a, b-1, c, d) = p;
                t_values(a, b-1, c, d) = stats.tstat;
            end
        end
    end
end
fprintf('done.\n')

% run permutations to construct null distribution 
fprintf('constructing null distribution: ')
for perm = 1:params.permutations
    % update
    if mod(perm, 100) == 0
        fprintf('. ')
    end

    % initiate output
    perm_p_values = NaN(size(p_values)); 
    perm_cluster_stats = [];

    % prepare shuffled dataset
    dataset_shuffled = dataset;
    for s = 1:size(dataset, 3)
        for a = 1:size(dataset, 1)
            for c = 1:size(dataset, 4)
                for b = 2:size(dataset, 2) 
                    % extract pre-patch and post-patch samples from the
                    % analysed time window
                    pre_patch = squeeze(dataset(a, 1, s, c, toi))';
                    post_patch = squeeze(dataset(a, b, s, c, toi))';

                    % combine data and shuffle labels
                    combined_samples = [pre_patch, post_patch];
                    shuffled_labels = randperm(length(combined_samples));
                
                    % reassign shuffled data
                    dataset_shuffled(a, 1, s, c, toi) = combined_samples(shuffled_labels(1:length(toi)));
                    dataset_shuffled(a, b, s, c, toi) = combined_samples(shuffled_labels(length(toi)+1:end));
                end
            end
        end
    end

    % recalculate t-tests on shuffled data
    for a = 1:size(dataset, 1)
        for b = 2:size(dataset, 2)
            for c = 1:size(dataset, 4)
                for d = 1:length(toi)
                    % extract pre-patch and post-patch values across subjects
                    pre_patch = squeeze(dataset_shuffled(a, 1, :, c, toi(d)));
                    post_patch = squeeze(dataset_shuffled(a, b, :, c, toi(d)));

                    % paired t-test on shuffled samples
                    [~, p, ~, stats] = ttest(post_patch, pre_patch);

                    % store results
                    perm_p_values(a, b-1, c, d) = p;
                    perm_t_values(a, b-1, c, d) = stats.tstat;
                end
            end
        end
    end

    % identify clusters in permuted data
    perm_cluster_stats = [];
    for a = 1:size(dataset, 1)
        for b = 1:size(dataset, 2)-1
            for c = 1:size(dataset, 4)
                % find significant time samples
                significant_samples = find(perm_p_values(a, b, c, :) < params.alpha);

                % compute cluster-level statistic = sum of t-values
                if ~isempty(significant_samples)
                    clusters = find_contiguous_clusters(significant_samples);           
                    for cluster = 1:length(clusters)
                        perm_cluster_stats = [perm_cluster_stats, sum(abs(perm_t_values(a, b, c, clusters{cluster})))];
                    end
                end
            end
        end
    end

    % store max cluster statistic from this permutation
    if ~isempty(perm_cluster_stats)
        max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        max_cluster_stats(perm) = 0;
    end
end
fprintf('done.\n')

% compute statistical threshold = 95th percentile of null distribution
threshold = prctile(max_cluster_stats, 95);

% identify significant clusters based on permutation results
fprintf('correcting p-values: ')
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        for c = 1:size(dataset, 4)
            % extract original clusters of p-values
            significant_samples = find(p_values(a, b, c, :) < params.alpha);
            if ~isempty(significant_samples)
                clusters = find_contiguous_clusters(significant_samples);
            
                % retain only clusters that pass the threshold
                if ~isempty(clusters)
                    for cluster = 1:length(clusters)
                        % compute overall 
                        cluster_t = sum(abs(t_values(a, b, c, clusters{cluster})));
                        if cluster_t > threshold
                            p_values_correct(a, b, c, clusters{cluster}) = 0.05; 
                        end
                    end
                end
            end
        end
    end
end
fprintf('done.\n')

% save data and measures to the output structures 
TEP_new_data.shuffled_time = dataset_shuffled;
TEP_new_measures.stats_time.toi = params.toi;
TEP_new_measures.stats_time.p_values = p_values;
TEP_new_measures.stats_time.t_values = t_values;
TEP_new_measures.stats_time.cluster_null_distribution = max_cluster_stats;
TEP_new_measures.stats_time.threshold = threshold;
row_counter = 1;
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        for c = 1:size(dataset, 4)
            % find significant samples
            significant_samples = find(p_values_correct(a, b, c, :) == 0.05);
            
            if ~isempty(significant_samples)                
                % identify contiguous clusters of significance
                clusters = find_contiguous_clusters(significant_samples);
                
                % store results
                for cluster = 1:length(clusters)
                    % convert sample indices to actual time values
                    significant_times = params.x((params.toi(1) - 0.001 + 1.5)*1000 + clusters{cluster});
                    
                    % encode
                    TEP_new_measures.stats_time.clusters(row_counter).condition = params.condition{a};
                    TEP_new_measures.stats_time.clusters(row_counter).timepoint = params.timepoint{b+1};
                    TEP_new_measures.stats_time.clusters(row_counter).electrode_n = c; 
                    TEP_new_measures.stats_time.clusters(row_counter).electrode_label = params.chanlocs(c).labels; 
                    TEP_new_measures.stats_time.clusters(row_counter).times{cluster} = significant_times;
                end

                % update row counter
                row_counter = row_counter + 1;
            end
        end
    end
end
fprintf('saving output:')
save(output_file, 'TEP_new_data', 'TEP_new_measures', '-append')
fprintf('done.\n')

% identify component TOIs
for c = 1:length(params.components)
    if strcmp(params.components{c}, 'N45')
        params.TOI(c, 1) = (TEP_new_measures.microstates.N45_P60.toi(1) + 0.5)/1000;
        params.TOI(c, 2) = (params.N45_offset - 0.5)/1000;
    elseif strcmp(params.components{c}, 'P60')
        params.TOI(c, 1) = (params.N45_offset + 0.5)/1000; 
        params.TOI(c, 2) = (TEP_new_measures.microstates.N45_P60.toi(2) - 0.5)/1000;
    else
        statement = sprintf('params.TOI(c, :) = TEP_new_measures.microstates.%s.toi;', params.components{c});  
        eval(statement)
        params.TOI(c, 1) = (params.TOI(c, 1) + 0.5)/1000;
        params.TOI(c, 2) = (params.TOI(c, 2) - 0.5)/1000;
    end
end

% plot significant clusters for each condition and post-patch timepoint
colors.palette = hsv(size(dataset, 4));
y.height = 0.8;
y.offset = 1:size(dataset, 4);
y.ticks = y.offset + y.height/2; 
y.labels = params.labels; 
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        % launch the figure
        fig = figure(figure_counter);
        screen_size = get(0, 'ScreenSize');
        set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 1.5])
        hold on

        % plot component TOIs
        for c = 1:length(params.components)
            rectangle('Position', [params.TOI(c, 1), y.ticks(1) - 1, params.TOI(c, 2) - params.TOI(c, 1), (y.ticks(end) + 1) - (y.ticks(1) - 1)], ...
                  'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none'); 
        end

        % select all electrodes containing significant clusters 
        idx = logical([]);
        for c = 1:length(TEP_new_measures.stats_time.clusters)
            if strcmp(TEP_new_measures.stats_time.clusters(c).condition, params.condition{a})
                if strcmp(TEP_new_measures.stats_time.clusters(c).timepoint, params.timepoint{b+1})
                    idx(c) = true;
                else
                    idx(c) = false;
                end
            else
                idx(c) = false;
            end
        end
        significant_electrodes = TEP_new_measures.stats_time.clusters(idx);

        % plot the clusters
        for e = 1:length(significant_electrodes)
            % assign y offset
            y.this = y.offset(significant_electrodes(e).electrode_n);

            % assign color
            colors.this = colors.palette(significant_electrodes(e).electrode_n, :);

            % loop through clusters
            for d = 1:length(significant_electrodes(e).times)                 
                % define rectangle position
                x.start = min(significant_electrodes(e).times{d});
                x.end = max(significant_electrodes(e).times{d});
                x.width = x.end - x.start;               
                    
                % plot rectangle
                rectangle('Position', [x.start, y.this, x.width, y.height], ...
                          'FaceColor', colors.this, 'EdgeColor', 'none'); 
            end
        end

        % adjust the figure
        title(sprintf('%s condition - %s', params.condition{a}, params.timepoint{b+1}));
        xlabel('time (s)');
        ylabel('electrodes');
        xlim([params.plot(1) params.plot(2)]);
        ylim([y.ticks(1) - 1, y.ticks(end) + 1])
        set(gca, 'FontSize', 14)
        
        % set y-ticks at center of rectangles
        yticks(y.ticks);
        yticklabels(y.labels);

        % axis visuals
        box off;
        ax = gca;
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';
        ax.TickDir = 'out'; 
        set(gca, 'Layer', 'Top')

        % plot stimulus
        y.limits = ylim;
        line([0, 0], y.limits, 'Color', [0.5020    0.5020    0.5020], 'LineWidth', 2.5, 'LineStyle', '--')

        % save figure and update counter
        saveas(fig, sprintf('%s\\figures\\change_%s_%s.svg', folder.output, params.condition{a}, params.timepoint{b+1}))
        figure_counter = figure_counter + 1;
    end
end

% plot average change across individual component time windows 
addpath(genpath([folder.toolbox '\EEGLAB']));
dataset = TEP_new_data.change;
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)
        % identify significant electrodes within this timepoint
        idx = logical([]);
        for d = 1:length(TEP_new_measures.stats_time.clusters)
            if strcmp(TEP_new_measures.stats_time.clusters(d).condition, params.condition{a})
                if strcmp(TEP_new_measures.stats_time.clusters(d).timepoint, params.timepoint{b+1})
                    idx(d) = true;
                else
                    idx(d) = false;
                end
            else
                idx(d) = false;
            end
        end
        significant_electrodes = TEP_new_measures.stats_time.clusters(idx);

        % cycle through components
        for c = 1:length(params.components)
            % identify x limits
            x.start = (params.TOI(c, 1) - 0.0005 - header.xstart)/header.xstep;
            x.end = (params.TOI(c, 2) + 0.0005 - header.xstart)/header.xstep;
    
            % calculate mean change
            for e = 1:size(dataset, 4)
                y.data(e) = mean(squeeze(dataset(a, b, :, e, x.start : x.end)), 'all');
            end

            % select significant electrodes
            chans2plot = [];
            if ~isempty(significant_electrodes)
                for d = 1:length(significant_electrodes)
                    for e = 1:length(significant_electrodes(d).times)
                        for f = significant_electrodes(d).times{e}
                            if f >= params.TOI(c, 1) && f <= params.TOI(c, 2)
                                chans2plot(end + 1) = significant_electrodes(d).electrode_n;
                            end
                        end
                    end
                end
            end
            
            % plot the topography
            fig = figure(figure_counter);
            topoplot(double(y.data'), params.chanlocs, 'style', 'map', 'shading', 'interp', 'whitebk', 'on', ...
                            'maplimits', params.topo_lim, 'electrodes', 'off');
            hold on
            title(sprintf('%s: %s - %s', params.components{c}, params.condition{a}, params.timepoint{b+1}))
            set(gca, 'FontSize', 14)

            % mark significant electrodes
            if ~isempty(chans2plot)
                for i = 1:length(chans2plot)
                    % extract electrode coordinates
                    electrode.theta = params.chanlocs(chans2plot(i)).theta - 90;   
                    electrode.radius = params.chanlocs(chans2plot(i)).radius;   

                    % convert polar to cartesian coordinates
                    electrode.x = electrode.radius .* cosd(electrode.theta);  
                    electrode.y = -electrode.radius .* sind(electrode.theta);  

                    % adjust for topoplot scaling 
                    electrode.scale_factor = 0.4 / 0.511;
                    electrode.x_scaled = electrode.x * electrode.scale_factor;
                    electrode.y_scaled = electrode.y * electrode.scale_factor;

                    hold on;
                    plot(electrode.x_scaled, electrode.y_scaled, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5);
                end
            end

            % % add colorbar
            % cb = colorbar('Location', 'southoutside');  
            % set(cb, 'Orientation', 'horizontal');  
            % ylabel(cb, 'amplitude (µV)', 'FontSize', 12);

            % save figure and update counter
            saveas(fig, sprintf('%s\\figures\\change_topo_%s_%s_%s.png', folder.output, params.components{c}, params.condition{a}, params.timepoint{b+1}))
            figure_counter = figure_counter + 1;
        end
    end
end

% clear and move on
clear a b c d e f p i header pre_patch post_patch stats p_values t_values toi ...
    combined_samples shuffled_labels dataset_shuffled p_values_correct perm_p_values perm_t_values perm_cluster_stats...
    threshold max_cluster_stats significant_samples clusters cluster cluster_t perm significant_times ...
    row_counter fig screen_size x y colors idx significant_electrodes ax statement data cb visual dataset electrode chans2plot 
fprintf('section 10 finished.\n')

%% 10b) statistics: change from baseline - spatiotemporal clusters
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.toi = [0.01 0.4];
params.permutations = 1000;
params.alpha = 0.05;
params.plot = [-0.025 0.4];
params.components = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
params.N45_offset = 58;
params.topo_lim = [-2, 2];
% -------------------------
fprintf('section 10: statistics - change from baseline\n')

% update
load(output_file, 'TEP_new_data', 'TEP_new_measures')
dataset = TEP_new_data.normalized_flipped;
cd(folder.output)
clear subject_idx

% load dataset information
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')
params.labels = {header.chanlocs.labels};
params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
params.chanlocs = header.chanlocs;

% identify analysed time window
toi = find(params.x >= params.toi(1) & params.x <= params.toi(2));

% initialize output
p_values = NaN(size(dataset, 1), size(dataset, 2) - 1, size(dataset, 4), length(toi));
p_values_correct = ones(size(p_values));
t_values = NaN(size(p_values));

% perform sample-wise paired t-test
fprintf('computing sample-wise paired t-test: ')
for a = 1:size(dataset, 1)
    for b = 2:size(dataset, 2)
        fprintf('. ')
        for c = 1:size(dataset, 4)
            for d = 1:length(toi)
                % extract pre-patch and post-patch values across subjects
                pre_patch = squeeze(dataset(a, 1, :, c, toi(d))); 
                post_patch = squeeze(dataset(a, b, :, c, toi(d))); 
                
                % perform paired t-test
                [~, p, ~, stats] = ttest(post_patch, pre_patch);
                
                % store results
                p_values(a, b-1, c, d) = p;
                t_values(a, b-1, c, d) = stats.tstat;
            end
        end
    end
end
fprintf('done.\n')

% run permutations to construct null distribution 
fprintf('constructing spatiotemporal null distribution: ')
for perm = 1:params.permutations
    % update
    if mod(perm, 100) == 0
        fprintf('. ')
    end

    % initiate output
    perm_p_values = NaN(size(p_values)); 
    max_cluster_stats(:, :, perm) = zeros(size(p_values, 1), size(p_values, 2));

    % prepare shuffled dataset
    dataset_shuffled = dataset;
    for s = 1:size(dataset, 3)
        for a = 1:size(dataset, 1)
            for c = 1:size(dataset, 4)
                for b = 2:size(dataset, 2) 
                    % extract pre-patch and post-patch samples from the
                    % analysed time window
                    pre_patch = squeeze(dataset(a, 1, s, c, toi))';
                    post_patch = squeeze(dataset(a, b, s, c, toi))';

                    % combine data and shuffle labels
                    combined_samples = [pre_patch, post_patch];
                    shuffled_labels = randperm(length(combined_samples));
                
                    % reassign shuffled data
                    dataset_shuffled(a, 1, s, c, toi) = combined_samples(shuffled_labels(1:length(toi)));
                    dataset_shuffled(a, b, s, c, toi) = combined_samples(shuffled_labels(length(toi)+1:end));
                end
            end
        end
    end

    % recalculate t-tests on shuffled data
    for a = 1:size(dataset, 1)
        for b = 2:size(dataset, 2)
            for c = 1:size(dataset, 4)
                for d = 1:length(toi)
                    % extract pre-patch and post-patch values across subjects
                    pre_patch = squeeze(dataset_shuffled(a, 1, :, c, toi(d)));
                    post_patch = squeeze(dataset_shuffled(a, b, :, c, toi(d)));

                    % paired t-test on shuffled samples
                    [~, p, ~, stats] = ttest(post_patch, pre_patch);

                    % store results
                    perm_p_values(a, b-1, c, d) = p;
                    perm_t_values(a, b-1, c, d) = stats.tstat;
                end
            end
        end
    end

    % identify spatiotemporal clusters in permuted data
    max_perm_cluster_stat = zeros(size(p_values, 1), size(p_values, 2));
    for a = 1:size(dataset, 1)
        for b = 1:size(dataset, 2)-1
            % find significant samples
            significant_matrix = squeeze(perm_p_values(a, b, :, :) < params.alpha);
            
            % find spatiotemporal clusters
            if any(significant_matrix(:))
                perm_clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);                
                for cluster = 1:length(perm_clusters)
                    cluster_indices = perm_clusters{cluster};
                    perm_cluster_stat = sum(abs(perm_t_values(a, b, cluster_indices(:,1), cluster_indices(:,2))), 'all');

                    % store max cluster statistic if eligible
                    if perm_cluster_stat > max_perm_cluster_stat(a, b)
                        max_perm_cluster_stat(a, b) = perm_cluster_stat;
                    end
                end
            end
        end
    end

    % store max cluster statistics from this permutation
    for a = 1:size(dataset, 1) 
        for b = 1:size(dataset, 2)-1 
            max_cluster_stats(a, b, perm) = max_perm_cluster_stat(a, b);
        end
    end
end
fprintf('done.\n')

% compute statistical thresholds = 95th percentile of null distribution
thresholds = prctile(max_cluster_stats, 95, 3);

% identify significant clusters based on permutation results
fprintf('identifying significant spatiotemporal clusters in real data: ');
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        % create binary matrix of significant time-electrode pairs
        significant_matrix = squeeze(p_values(a, b, :, :) < params.alpha);

        % retain only clusters that pass the threshold
        if any(significant_matrix(:))
            clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);
            if ~isempty(clusters)
                for cluster = 1:length(clusters)
                    % compute overall cluster t-statistic
                    cluster_indices = clusters{cluster};
                    cluster_t = sum(abs(t_values(a, b, cluster_indices(:,1), cluster_indices(:,2))));
                
                    % store cluster if significant 
                    if cluster_t > thresholds(a, b)
                        for i = 1:size(cluster_indices, 1)
                            p_values_correct(a, b, cluster_indices(i,1), cluster_indices(i,2)) = 0.05;
                        end
                    end
                end
            end
        end
    end
end
fprintf('done.\n')

% save data and measures to the output structures 
TEP_new_data.shuffled_time_ST = dataset_shuffled;
TEP_new_measures.stats_time_ST.toi = params.toi;
TEP_new_measures.stats_time_ST.p_values = p_values;
TEP_new_measures.stats_time_ST.t_values = t_values;
TEP_new_measures.stats_time_ST.cluster_null_distributions = max_cluster_stats;
TEP_new_measures.stats_time_ST.thresholds = thresholds;
row_counter = 1;
row_counter_ST = 1;
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        % find significant samples per electrode
        for c = 1:size(dataset, 4)            
            significant_samples = find(p_values_correct(a, b, c, :) == 0.05);            
            if ~isempty(significant_samples)                
                % identify contiguous clusters of significance
                clusters = find_contiguous_clusters(significant_samples);
                
                % store results
                for cluster = 1:length(clusters)
                    % convert sample indices to actual time values
                    significant_times = params.x((params.toi(1) - 0.001 + 1.5)*1000 + clusters{cluster});
                    
                    % encode
                    TEP_new_measures.stats_time_ST.clusters(row_counter).condition = params.condition{a};
                    TEP_new_measures.stats_time_ST.clusters(row_counter).timepoint = params.timepoint{b+1};
                    TEP_new_measures.stats_time_ST.clusters(row_counter).electrode_n = c; 
                    TEP_new_measures.stats_time_ST.clusters(row_counter).electrode_label = params.chanlocs(c).labels; 
                    TEP_new_measures.stats_time_ST.clusters(row_counter).times{cluster} = significant_times;
                end

                % update row counter
                row_counter = row_counter + 1;
            end
        end

        % find significant spatiotemporal clusters
        significant_matrix = squeeze(p_values_correct(a, b, :, :) == 0.05);
        if any(significant_matrix(:))
            clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);
            if ~isempty(clusters)
                % encode
                TEP_new_measures.stats_time_ST.clusters_ST(row_counter_ST).condition = params.condition{a};
                TEP_new_measures.stats_time_ST.clusters_ST(row_counter_ST).timepoint = params.timepoint{b+1};
                TEP_new_measures.stats_time_ST.clusters_ST(row_counter_ST).STclusters = clusters;
            end

            % update row counter
            row_counter_ST = row_counter_ST + 1;
        end
    end
end
fprintf('saving output:')
save(output_file, 'TEP_new_data', 'TEP_new_measures', '-append')
fprintf('done.\n')

% plot null distributions & thresholds
figure('units','normalized','outerposition',[0 0 1 1]);
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        % select data
        visual.values = squeeze(TEP_new_measures.stats_time_ST.cluster_null_distributions(a, b, :))';  
        visual.threshold = TEP_new_measures.stats_time_ST.thresholds(a, b);
        visual.title = sprintf('%s condition, %s', params.condition{a}, params.timepoint{b + 1});

        % choose colour 
        if a == 1
            visual.color = [0.9216    0.1490    0.1490];
        elseif a == 2
            visual.color = [0    0.4471    0.7412];
        end

        % plot 
        subplot(size(dataset, 1), size(dataset, 2)-1, (a-1)*(size(dataset, 2)-1) + b)
        hold on
        plot_distribution(visual)
    end
end
saveas(fig, sprintf('%s\\figures\\stats_time_ST_distribution.png', folder.output))
figure_counter = figure_counter + 1;

% identify component TOIs
for c = 1:length(params.components)
    if strcmp(params.components{c}, 'N45')
        params.TOI(c, 1) = (TEP_new_measures.microstates.N45_P60.toi(1) + 0.5)/1000;
        params.TOI(c, 2) = (params.N45_offset - 0.5)/1000;
    elseif strcmp(params.components{c}, 'P60')
        params.TOI(c, 1) = (params.N45_offset + 0.5)/1000; 
        params.TOI(c, 2) = (TEP_new_measures.microstates.N45_P60.toi(2) - 0.5)/1000;
    else
        statement = sprintf('params.TOI(c, :) = TEP_new_measures.microstates.%s.toi;', params.components{c});  
        eval(statement)
        params.TOI(c, 1) = (params.TOI(c, 1) + 0.5)/1000;
        params.TOI(c, 2) = (params.TOI(c, 2) - 0.5)/1000;
    end
end

% plot significant clusters for each condition and post-patch timepoint
colors.palette = hsv(size(dataset, 4));
y.height = 0.8;
y.offset = 1:size(dataset, 4);
y.ticks = y.offset + y.height/2; 
y.labels = params.labels; 
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)-1
        % launch the figure
        fig = figure(figure_counter);
        screen_size = get(0, 'ScreenSize');
        set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 1.5])
        hold on

        % plot component TOIs
        for c = 1:length(params.components)
            rectangle('Position', [params.TOI(c, 1), y.ticks(1) - 1, params.TOI(c, 2) - params.TOI(c, 1), (y.ticks(end) + 1) - (y.ticks(1) - 1)], ...
                  'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none'); 
        end

        % select all electrodes containing significant clusters 
        idx = logical([]);
        for c = 1:length(TEP_new_measures.stats_time_ST.clusters)
            if strcmp(TEP_new_measures.stats_time_ST.clusters(c).condition, params.condition{a})
                if strcmp(TEP_new_measures.stats_time_ST.clusters(c).timepoint, params.timepoint{b+1})
                    idx(c) = true;
                else
                    idx(c) = false;
                end
            else
                idx(c) = false;
            end
        end
        significant_electrodes = TEP_new_measures.stats_time_ST.clusters(idx);

        % plot the clusters
        for e = 1:length(significant_electrodes)
            % assign y offset
            y.this = y.offset(significant_electrodes(e).electrode_n);

            % assign color
            colors.this = colors.palette(significant_electrodes(e).electrode_n, :);

            % loop through clusters
            for d = 1:length(significant_electrodes(e).times)                 
                % define rectangle position
                x.start = min(significant_electrodes(e).times{d});
                x.end = max(significant_electrodes(e).times{d});
                x.width = x.end - x.start;               
                    
                % plot rectangle
                rectangle('Position', [x.start, y.this, x.width, y.height], ...
                          'FaceColor', colors.this, 'EdgeColor', 'none'); 
            end
        end

        % adjust the figure
        title(sprintf('%s condition - %s', params.condition{a}, params.timepoint{b+1}));
        xlabel('time (s)');
        ylabel('electrodes');
        xlim([params.plot(1) params.plot(2)]);
        ylim([y.ticks(1) - 1, y.ticks(end) + 1])
        set(gca, 'FontSize', 14)
        
        % set y-ticks at center of rectangles
        yticks(y.ticks);
        yticklabels(y.labels);

        % axis visuals
        box off;
        ax = gca;
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';
        ax.TickDir = 'out'; 
        set(gca, 'Layer', 'Top')

        % plot stimulus
        y.limits = ylim;
        line([0, 0], y.limits, 'Color', [0.5020    0.5020    0.5020], 'LineWidth', 2.5, 'LineStyle', '--')

        % save figure and update counter
        saveas(fig, sprintf('%s\\figures\\change_%s_%s.svg', folder.output, params.condition{a}, params.timepoint{b+1}))
        figure_counter = figure_counter + 1;
    end
end

% plot average change across individual component time windows 
addpath(genpath([folder.toolbox '\EEGLAB']));
dataset = TEP_new_data.change;
for a = 1:size(dataset, 1)
    for b = 1:size(dataset, 2)
        % identify significant electrodes within this timepoint
        idx = logical([]);
        for d = 1:length(TEP_new_measures.stats_time_ST.clusters)
            if strcmp(TEP_new_measures.stats_time_ST.clusters(d).condition, params.condition{a})
                if strcmp(TEP_new_measures.stats_time_ST.clusters(d).timepoint, params.timepoint{b+1})
                    idx(d) = true;
                else
                    idx(d) = false;
                end
            else
                idx(d) = false;
            end
        end
        significant_electrodes = TEP_new_measures.stats_time_ST.clusters(idx);

        % cycle through components
        for c = 1:length(params.components)
            % identify x limits
            x.start = (params.TOI(c, 1) - 0.0005 - header.xstart)/header.xstep;
            x.end = (params.TOI(c, 2) + 0.0005 - header.xstart)/header.xstep;
    
            % calculate mean change
            for e = 1:size(dataset, 4)
                y.data(e) = mean(squeeze(dataset(a, b, :, e, x.start : x.end)), 'all');
            end

            % select significant electrodes
            chans2plot = [];
            if ~isempty(significant_electrodes)
                for d = 1:length(significant_electrodes)
                    for e = 1:length(significant_electrodes(d).times)
                        for f = significant_electrodes(d).times{e}
                            if f >= params.TOI(c, 1) && f <= params.TOI(c, 2)
                                chans2plot(end + 1) = significant_electrodes(d).electrode_n;
                            end
                        end
                    end
                end
            end
            
            % plot the topography
            fig = figure(figure_counter);
            topoplot(double(y.data'), params.chanlocs, 'style', 'map', 'shading', 'interp', 'whitebk', 'on', ...
                            'maplimits', params.topo_lim, 'electrodes', 'off');
            hold on
            title(sprintf('%s: %s - %s', params.components{c}, params.condition{a}, params.timepoint{b+1}))
            set(gca, 'FontSize', 14)

            % mark significant electrodes
            if ~isempty(chans2plot)
                for i = 1:length(chans2plot)
                    % extract electrode coordinates
                    electrode.theta = params.chanlocs(chans2plot(i)).theta - 90;   
                    electrode.radius = params.chanlocs(chans2plot(i)).radius;   

                    % convert polar to cartesian coordinates
                    electrode.x = electrode.radius .* cosd(electrode.theta);  
                    electrode.y = -electrode.radius .* sind(electrode.theta);  

                    % adjust for topoplot scaling 
                    electrode.scale_factor = 0.4 / 0.511;
                    electrode.x_scaled = electrode.x * electrode.scale_factor;
                    electrode.y_scaled = electrode.y * electrode.scale_factor;

                    hold on;
                    plot(electrode.x_scaled, electrode.y_scaled, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5);
                end
            end

            % % add colorbar
            % cb = colorbar('Location', 'southoutside');  
            % set(cb, 'Orientation', 'horizontal');  
            % ylabel(cb, 'amplitude (µV)', 'FontSize', 12);

            % save figure and update counter
            saveas(fig, sprintf('%s\\figures\\change_topo_%s_%s_%s.png', folder.output, params.components{c}, params.condition{a}, params.timepoint{b+1}))
            figure_counter = figure_counter + 1;
        end
    end
end

% clear and move on
clear a b c d e i j f p s header pre_patch post_patch stats p_values t_values toi neighbors ...
    combined_samples shuffled_labels dataset_shuffled p_values_correct perm_p_values perm_t_values perm_clusters perm_cluster_stat significant_matrix ...
    max_perm_cluster_stat max_cluster_stats thresholds significant_samples clusters cluster cluster_t perm significant_times cluster_indices ...
    row_counter row_counter_ST fig screen_size x y colors idx significant_electrodes ax statement data cb visual dataset electrode chans2plot 
fprintf('section 10 finished.\n')

%% 11a) statistics: differences across conditions - separately per electrode
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.toi = [0.01 0.4];
params.max_distance = 0.7;
params.permutations = 1000;
params.alpha = 0.05;
params.plot = [-0.025 0.4];
params.components = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
params.N45_offset = 58;
params.topo_lim = [-2, 2];
% -------------------------
fprintf('section 11: statistics - differences across conditions and time\n')

% update
load(output_file, 'TEP_new_data', 'TEP_new_measures')
dataset = TEP_new_data.change;
cd(folder.output)
clear subject_idx

% load dataset information
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')
params.labels = {header.chanlocs.labels};
params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
params.chanlocs = header.chanlocs;

% identify analysed time window
toi = find(params.x >= params.toi(1) & params.x <= params.toi(2));

% initialize output
p_values = NaN(size(dataset, 2), size(dataset, 4), length(toi));
p_values_correct = ones(size(p_values));
t_values = NaN(size(p_values));

% identify neighboring electrodes
neighbors = find_electrode_neighbors(header.chanlocs, params.max_distance);

% perform sample-wise paired t-test
fprintf('computing sample-wise paired t-test: ')
for b = 1:size(dataset, 2)
    fprintf('. ')
    for c = 1:size(dataset, 4)
        for d = 1:length(toi)
            % extract data from both conditions across subjects
            data_pain = squeeze(dataset(1, b, :, c, toi(d)))'; 
            data_control = squeeze(dataset(2, b, :, c, toi(d)))'; 
            
            % perform paired t-test
            [~, p, ~, stats] = ttest(data_pain, data_control);
            
            % store results
            p_values(b, c, d) = p;
            t_values(b, c, d) = stats.tstat;
        end
    end
end
fprintf('done.\n')

% run permutations to construct null distribution 
fprintf('constructing null distribution: ')
for perm = 1:params.permutations
    % update
    if mod(perm, 100) == 0
        fprintf('. ')
    end

    % initiate output
    perm_p_values = NaN(size(p_values)); 
    perm_cluster_stats = [];

    % prepare shuffled dataset
    dataset_shuffled = dataset;
    for s = 1:size(dataset, 3)
        for c = 1:size(dataset, 4)
            for b = 1:size(dataset, 2) 
                % extract data from both conditions across the analysed toi
                data_pain = squeeze(dataset(1, b, s, c, toi))'; 
                data_control = squeeze(dataset(2, b, s, c, toi))'; 

                % combine data and shuffle labels
                combined_samples = [data_control, data_pain];
                shuffled_labels = randperm(length(combined_samples));
            
                % reassign shuffled data
                dataset_shuffled(1, b, s, c, toi) = combined_samples(shuffled_labels(1:length(toi)));
                dataset_shuffled(2, b, s, c, toi) = combined_samples(shuffled_labels(length(toi)+1:end));
            end
        end
    end

    % recalculate t-tests on shuffled data
    for b = 1:size(dataset, 2)
        for c = 1:size(dataset, 4)
            for d = 1:length(toi)
                % extract pre-patch and post-patch values across subjects
                data_pain = squeeze(dataset_shuffled(1, b, :, c, toi(d)))'; 
                data_control = squeeze(dataset_shuffled(2, b, :, c, toi(d)))'; 

                % paired t-test on shuffled samples
                [~, p, ~, stats] = ttest(data_pain, data_control);

                % store results
                perm_p_values(b, c, d) = p;
                perm_t_values(b, c, d) = stats.tstat;
            end
        end
    end

    % identify clusters in permuted data
    perm_cluster_stats = [];
    for b = 1:size(dataset, 2)
        for c = 1:size(dataset, 4)
            % find significant time samples
            significant_samples = find(perm_p_values(b, c, :) < params.alpha);

            % compute cluster-level statistic = sum of t-values
            if ~isempty(significant_samples)
                clusters = find_contiguous_clusters(significant_samples);           
                for cluster = 1:length(clusters)
                    perm_cluster_stats = [perm_cluster_stats, sum(abs(perm_t_values(b, c, clusters{cluster})))];
                end
            end
        end
    end

    % store max cluster statistic from this permutation
    if ~isempty(perm_cluster_stats)
        max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        max_cluster_stats(perm) = 0;
    end
end
fprintf('done.\n')

% compute statistical threshold = 95th percentile of null distribution
threshold = prctile(max_cluster_stats, 95);

% identify significant clusters based on permutation results
fprintf('correcting p-values: ')
for b = 1:size(dataset, 2)
    for c = 1:size(dataset, 4)
        % extract original clusters of p-values
        significant_samples = find(p_values(b, c, :) < params.alpha);
        if ~isempty(significant_samples)
            clusters = find_contiguous_clusters(significant_samples);
        
            % retain only clusters that pass the threshold
            if ~isempty(clusters)
                for cluster = 1:length(clusters)
                    % compute overall 
                    cluster_t = sum(abs(t_values(b, c, clusters{cluster})));
                    if cluster_t > threshold
                        p_values_correct(b, c, clusters{cluster}) = 0.05; 
                    end
                end
            end
        end
    end
end
fprintf('done.\n')

% save data and measures to the output structures 
TEP_new_data.shuffled_condition = dataset_shuffled;
TEP_new_measures.stats_condition.toi = params.toi;
TEP_new_measures.stats_condition.p_values = p_values;
TEP_new_measures.stats_condition.t_values = t_values;
TEP_new_measures.stats_condition.cluster_null_distribution = max_cluster_stats;
TEP_new_measures.stats_condition.threshold = threshold;
row_counter = 1;
for b = 1:size(dataset, 2)
    for c = 1:size(dataset, 4)
        % find significant samples
        significant_samples = find(p_values_correct(b, c, :) == 0.05);
        
        if ~isempty(significant_samples)                
            % identify contiguous clusters of significance
            clusters = find_contiguous_clusters(significant_samples);
            
            % store results
            for cluster = 1:length(clusters)
                % convert sample indices to actual time values
                significant_times = params.x((params.toi(1) - 0.001 + 1.5)*1000 + clusters{cluster});
                
                % encode
                TEP_new_measures.stats_condition.clusters(row_counter).timepoint = params.timepoint{b+1};
                TEP_new_measures.stats_condition.clusters(row_counter).electrode_n = c; 
                TEP_new_measures.stats_condition.clusters(row_counter).electrode_label = params.chanlocs(c).labels; 
                TEP_new_measures.stats_condition.clusters(row_counter).times{cluster} = significant_times;
            end

            % update row counter
            row_counter = row_counter + 1;
        end
    end
end
fprintf('saving output:')
save(output_file, 'TEP_new_data', 'TEP_new_measures', '-append')
fprintf('done.\n')

% identify component TOIs
for c = 1:length(params.components)
    if strcmp(params.components{c}, 'N45')
        params.TOI(c, 1) = (TEP_new_measures.microstates.N45_P60.toi(1) + 0.5)/1000;
        params.TOI(c, 2) = (params.N45_offset - 0.5)/1000;
    elseif strcmp(params.components{c}, 'P60')
        params.TOI(c, 1) = (params.N45_offset + 0.5)/1000; 
        params.TOI(c, 2) = (TEP_new_measures.microstates.N45_P60.toi(2) - 0.5)/1000;
    else
        statement = sprintf('params.TOI(c, :) = TEP_new_measures.microstates.%s.toi;', params.components{c});  
        eval(statement)
        params.TOI(c, 1) = (params.TOI(c, 1) + 0.5)/1000;
        params.TOI(c, 2) = (params.TOI(c, 2) - 0.5)/1000;
    end
end

% plot significant clusters for each condition and post-patch timepoint
colors.palette = hsv(size(dataset, 4));
y.height = 0.8;
y.offset = 1:size(dataset, 4);
y.ticks = y.offset + y.height/2; 
y.labels = params.labels; 
for b = 1:size(dataset, 2)
    % launch the figure
    fig = figure(figure_counter);
    screen_size = get(0, 'ScreenSize');
    set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 1.5])
    hold on

    % plot component TOIs
    for c = 1:length(params.components)
        rectangle('Position', [params.TOI(c, 1), y.ticks(1) - 1, params.TOI(c, 2) - params.TOI(c, 1), (y.ticks(end) + 1) - (y.ticks(1) - 1)], ...
              'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none'); 
    end

    % select all electrodes containing significant clusters 
    idx = logical([]);
    for c = 1:length(TEP_new_measures.stats_condition.clusters)
        if strcmp(TEP_new_measures.stats_condition.clusters(c).timepoint, params.timepoint{b+1})
            idx(c) = true;
        else
            idx(c) = false;
        end
    end
    significant_electrodes = TEP_new_measures.stats_condition.clusters(idx);

    % plot the clusters
    for e = 1:length(significant_electrodes)
        % assign y offset
        y.this = y.offset(significant_electrodes(e).electrode_n);

        % assign color
        colors.this = colors.palette(significant_electrodes(e).electrode_n, :);

        % loop through clusters
        for d = 1:length(significant_electrodes(e).times)                 
            % define rectangle position
            x.start = min(significant_electrodes(e).times{d});
            x.end = max(significant_electrodes(e).times{d});
            x.width = x.end - x.start;               
                
            % plot rectangle
            rectangle('Position', [x.start, y.this, x.width, y.height], ...
                      'FaceColor', colors.this, 'EdgeColor', 'none'); 
        end
    end

    % adjust the figure
    title(sprintf('interaction - %s', params.timepoint{b+1}));
    xlabel('time (s)');
    ylabel('electrodes');
    xlim([params.plot(1) params.plot(2)]);
    ylim([y.ticks(1) - 1, y.ticks(end) + 1])
    set(gca, 'FontSize', 14)
    
    % set y-ticks at center of rectangles
    yticks(y.ticks);
    yticklabels(y.labels);

    % axis visuals
    box off;
    ax = gca;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.TickDir = 'out'; 
    set(gca, 'Layer', 'Top')

    % plot stimulus
    y.limits = ylim;
    line([0, 0], y.limits, 'Color', [0.5020    0.5020    0.5020], 'LineWidth', 2.5, 'LineStyle', '--')

    % save figure and update counter
    saveas(fig, sprintf('%s\\figures\\change_interaction_%s.svg', folder.output, params.timepoint{b+1}))
    figure_counter = figure_counter + 1;
end

% clear and move on
clear a b c d s p header dataset dataset_shuffled toi p_values p_values_correct t_values data_pain data_control ....
        perm perm_p_values perm_cluster_stats combined_samples shuffled_labels cluster significant_samples max_cluster_stats ...
        row_counter fig screen_size x y colors idx significant_electrodes ax statement data visual dataset 
fprintf('section 11 finished.\n')

%% 11b) statistics: differences across conditions - spatiotemporal clusters
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.toi = [0.01 0.4];
params.max_distance = 0.7;
params.permutations = 1000;
params.alpha = 0.05;
params.plot = [-0.025 0.4];
params.components = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
params.N45_offset = 58;
params.topo_lim = [-2, 2];
% -------------------------
fprintf('section 11: statistics - differences across conditions and time\n')

% update
load(output_file, 'TEP_new_data', 'TEP_new_measures')
dataset = TEP_new_data.change;
cd(folder.output)
clear subject_idx

% load dataset information
load(sprintf('%s\\%s %s %s %s %s .lw6', folder.processed, params.prefix, study, TEP_new(1).ID, params.condition{1}, params.timepoint{1}), '-mat')
params.labels = {header.chanlocs.labels};
params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
params.chanlocs = header.chanlocs;

% identify analysed time window
toi = find(params.x >= params.toi(1) & params.x <= params.toi(2));

% initialize output
p_values = NaN(size(dataset, 2), size(dataset, 4), length(toi));
p_values_correct = ones(size(p_values));
t_values = NaN(size(p_values));

% identify neighboring electrodes
neighbors = find_electrode_neighbors(header.chanlocs, params.max_distance);

% perform sample-wise paired t-test
fprintf('computing sample-wise paired t-test: ')
for b = 1:size(dataset, 2)
    fprintf('. ')
    for c = 1:size(dataset, 4)
        for d = 1:length(toi)
            % extract data from both conditions across subjects
            data_pain = squeeze(dataset(1, b, :, c, toi(d)))'; 
            data_control = squeeze(dataset(2, b, :, c, toi(d)))'; 
            
            % perform paired t-test
            [~, p, ~, stats] = ttest(data_pain, data_control);
            
            % store results
            p_values(b, c, d) = p;
            t_values(b, c, d) = stats.tstat;
        end
    end
end
fprintf('done.\n')

% run permutations to construct null distribution 
fprintf('constructing spatiotemporal null distribution: ')
for perm = 1:params.permutations
    % update
    if mod(perm, 100) == 0
        fprintf('. ')
    end

    % initiate output
    perm_p_values = NaN(size(p_values)); 
    max_cluster_stats(:, perm) = zeros(size(p_values, 1), 1);

    % prepare shuffled dataset
    dataset_shuffled = dataset;
    for s = 1:size(dataset, 3)
        for c = 1:size(dataset, 4)
            for b = 1:size(dataset, 2) 
                % extract data from both conditions across the analysed toi
                data_pain = squeeze(dataset(1, b, s, c, toi))'; 
                data_control = squeeze(dataset(2, b, s, c, toi))'; 

                % combine data and shuffle labels
                combined_samples = [data_control, data_pain];
                shuffled_labels = randperm(length(combined_samples));
            
                % reassign shuffled data
                dataset_shuffled(1, b, s, c, toi) = combined_samples(shuffled_labels(1:length(toi)));
                dataset_shuffled(2, b, s, c, toi) = combined_samples(shuffled_labels(length(toi)+1:end));
            end
        end
    end

    % recalculate t-tests on shuffled data
    for b = 1:size(dataset, 2)
        for c = 1:size(dataset, 4)
            for d = 1:length(toi)
                % extract pre-patch and post-patch values across subjects
                data_pain = squeeze(dataset_shuffled(1, b, :, c, toi(d)))'; 
                data_control = squeeze(dataset_shuffled(2, b, :, c, toi(d)))'; 

                % paired t-test on shuffled samples
                [~, p, ~, stats] = ttest(data_pain, data_control);

                % store results
                perm_p_values(b, c, d) = p;
                perm_t_values(b, c, d) = stats.tstat;
            end
        end
    end

    % identify spatiotemporal clusters in permuted data
    max_perm_cluster_stat = zeros(1, size(p_values, 1));
    for b = 1:size(dataset, 2)
        % create binary matrix of significant time-electrode pairs
        significant_matrix = squeeze(perm_p_values(b, :, :) < params.alpha);

        % find spatiotemporal clusters
        if any(significant_matrix(:))
            perm_clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);                
            for cluster = 1:length(perm_clusters)
                cluster_indices = perm_clusters{cluster};
                perm_cluster_stat = sum(abs(perm_t_values(b, cluster_indices(:,1), cluster_indices(:,2))), 'all');

                % store max cluster statistic if eligible
                if perm_cluster_stat > max_perm_cluster_stat(b)
                    max_perm_cluster_stat(b) = perm_cluster_stat;
                end
            end
        end
    end

    % store max cluster statistic from this permutation
    for b = 1:size(dataset, 2)
        max_cluster_stats(b, perm) = max_perm_cluster_stat(b);
    end
end
fprintf('done.\n')

% compute statistical thresholds = 95th percentile of null distribution
thresholds = prctile(max_cluster_stats, 95, 2);

% identify significant clusters based on permutation results
fprintf('identifying significant spatiotemporal clusters in real data: ')
for b = 1:size(dataset, 2)
    % create binary matrix of significant time-electrode pairs
    significant_matrix = squeeze(p_values(b, :, :) < params.alpha);

    % retain only clusters that pass the threshold
    if any(significant_matrix(:))
        clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);
        if ~isempty(clusters)
            for cluster = 1:length(clusters)
                % compute overall cluster t-statistic
                cluster_indices = clusters{cluster};
                cluster_t = sum(abs(t_values(b, cluster_indices(:,1), cluster_indices(:,2))));
            
                % store cluster if significant 
                if cluster_t > thresholds(b)
                    for i = 1:size(cluster_indices, 1)
                        p_values_correct(b, cluster_indices(i,1), cluster_indices(i,2)) = 0.05;
                    end
                end
            end
        end
    end
end
fprintf('done.\n')

% save data and measures to the output structures 
TEP_new_data.shuffled_condition_ST = dataset_shuffled;
TEP_new_measures.stats_condition_ST.toi = params.toi;
TEP_new_measures.stats_condition_ST.p_values = p_values;
TEP_new_measures.stats_condition_ST.t_values = t_values;
TEP_new_measures.stats_condition_ST.cluster_null_distribution = max_cluster_stats;
TEP_new_measures.stats_condition_ST.threshold = thresholds;
row_counter = 1;
row_counter_ST = 1;
for b = 1:size(dataset, 2)
    % find significant samples per electrode
    for c = 1:size(dataset, 4)
        significant_samples = find(p_values_correct(b, c, :) == 0.05);        
        if ~isempty(significant_samples)                
            % identify contiguous clusters of significance
            clusters = find_contiguous_clusters(significant_samples);
            
            % store results
            for cluster = 1:length(clusters)
                % convert sample indices to actual time values
                significant_times = params.x((params.toi(1) - 0.001 + 1.5)*1000 + clusters{cluster});
                
                % encode
                TEP_new_measures.stats_condition_ST.clusters(row_counter).timepoint = params.timepoint{b+1};
                TEP_new_measures.stats_condition_ST.clusters(row_counter).electrode_n = c; 
                TEP_new_measures.stats_condition_ST.clusters(row_counter).electrode_label = params.chanlocs(c).labels; 
                TEP_new_measures.stats_condition_ST.clusters(row_counter).times{cluster} = significant_times;
            end

            % update row counter
            row_counter = row_counter + 1;
        end
    end

    % find significant spatiotemporal clusters
    significant_matrix = squeeze(p_values_correct(b, :, :) == 0.05);
    if any(significant_matrix(:))
        clusters = find_spatiotemporal_clusters(significant_matrix, neighbors);
        if ~isempty(clusters)
            % encode
            TEP_new_measures.stats_condition_ST.clusters_ST(row_counter_ST).timepoint = params.timepoint{b+1};
            TEP_new_measures.stats_condition_ST.clusters_ST(row_counter_ST).STclusters = clusters;
        end

        % update row counter
        row_counter_ST = row_counter_ST + 1;
    end
end
fprintf('saving output:')
save(output_file, 'TEP_new_data', 'TEP_new_measures', '-append')
fprintf('done.\n')

% plot null distributions & thresholds
figure('units','normalized','outerposition',[0 0 1 1]);
for b = 1:size(dataset, 2)
    % select data
    visual.values = squeeze(TEP_new_measures.stats_condition_ST.cluster_null_distributions(b, :))';  
    visual.threshold = TEP_new_measures.stats_condition_ST.thresholds(b);
    visual.title = sprintf('%s', params.timepoint{b + 1});
    visual.color = [1.0000    0.0745    0.6510];

    % plot 
    subplot(2, size(dataset, 2)/2, b)
    hold on
    plot_distribution(visual)
end
saveas(gcf, sprintf('%s\\figures\\stats_condition_ST_distribution.png', folder.output))
figure_counter = figure_counter + 1;

% identify component TOIs
for c = 1:length(params.components)
    if strcmp(params.components{c}, 'N45')
        params.TOI(c, 1) = (TEP_new_measures.microstates.N45_P60.toi(1) + 0.5)/1000;
        params.TOI(c, 2) = (params.N45_offset - 0.5)/1000;
    elseif strcmp(params.components{c}, 'P60')
        params.TOI(c, 1) = (params.N45_offset + 0.5)/1000; 
        params.TOI(c, 2) = (TEP_new_measures.microstates.N45_P60.toi(2) - 0.5)/1000;
    else
        statement = sprintf('params.TOI(c, :) = TEP_new_measures.microstates.%s.toi;', params.components{c});  
        eval(statement)
        params.TOI(c, 1) = (params.TOI(c, 1) + 0.5)/1000;
        params.TOI(c, 2) = (params.TOI(c, 2) - 0.5)/1000;
    end
end

% plot significant clusters for each condition and post-patch timepoint
colors.palette = hsv(size(dataset, 4));
y.height = 0.8;
y.offset = 1:size(dataset, 4);
y.ticks = y.offset + y.height/2; 
y.labels = params.labels; 
for b = 1:size(dataset, 2)
    % launch the figure
    fig = figure(figure_counter);
    screen_size = get(0, 'ScreenSize');
    set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3)/2, screen_size(4) / 1.5])
    hold on

    % plot component TOIs
    for c = 1:length(params.components)
        rectangle('Position', [params.TOI(c, 1), y.ticks(1) - 1, params.TOI(c, 2) - params.TOI(c, 1), (y.ticks(end) + 1) - (y.ticks(1) - 1)], ...
              'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none'); 
    end

    % select all electrodes containing significant clusters 
    idx = logical([]);
    for c = 1:length(TEP_new_measures.stats_condition_ST.clusters)
        if strcmp(TEP_new_measures.stats_condition_ST.clusters(c).timepoint, params.timepoint{b+1})
            idx(c) = true;
        else
            idx(c) = false;
        end
    end
    significant_electrodes = TEP_new_measures.stats_condition_ST.clusters(idx);

    % plot the clusters
    for e = 1:length(significant_electrodes)
        % assign y offset
        y.this = y.offset(significant_electrodes(e).electrode_n);

        % assign color
        colors.this = colors.palette(significant_electrodes(e).electrode_n, :);

        % loop through clusters
        for d = 1:length(significant_electrodes(e).times)                 
            % define rectangle position
            x.start = min(significant_electrodes(e).times{d});
            x.end = max(significant_electrodes(e).times{d});
            x.width = x.end - x.start;               
                
            % plot rectangle
            rectangle('Position', [x.start, y.this, x.width, y.height], ...
                      'FaceColor', colors.this, 'EdgeColor', 'none'); 
        end
    end

    % adjust the figure
    title(sprintf('interaction - %s', params.timepoint{b+1}));
    xlabel('time (s)');
    ylabel('electrodes');
    xlim([params.plot(1) params.plot(2)]);
    ylim([y.ticks(1) - 1, y.ticks(end) + 1])
    set(gca, 'FontSize', 14)
    
    % set y-ticks at center of rectangles
    yticks(y.ticks);
    yticklabels(y.labels);

    % axis visuals
    box off;
    ax = gca;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.TickDir = 'out'; 
    set(gca, 'Layer', 'Top')

    % plot stimulus
    y.limits = ylim;
    line([0, 0], y.limits, 'Color', [0.5020    0.5020    0.5020], 'LineWidth', 2.5, 'LineStyle', '--')

    % save figure and update counter
    saveas(fig, sprintf('%s\\figures\\change_interaction_%s.svg', folder.output, params.timepoint{b+1}))
    figure_counter = figure_counter + 1;
end

% plot average change across individual component time windows 
addpath(genpath([folder.toolbox '\EEGLAB']));
for b = 1:size(dataset, 2)
    % identify significant electrodes within this timepoint
    idx = logical([]);
    for d = 1:length(TEP_new_measures.stats_condition_ST.clusters)
        if strcmp(TEP_new_measures.stats_condition_ST.clusters(d).timepoint, params.timepoint{b+1})
            idx(d) = true;
        else
            idx(d) = false;
        end
    end
    significant_electrodes = TEP_new_measures.stats_condition_ST.clusters(idx);

    % cycle through components
    for c = 1:length(params.components)
        % identify x limits
        x.start = (params.TOI(c, 1) - 0.0005 - header.xstart)/header.xstep;
        x.end = (params.TOI(c, 2) + 0.0005 - header.xstart)/header.xstep;

        % calculate mean change
        y.data = [];
        for s = 1:size(dataset, 3)
            for e = 1:size(dataset, 4)
                y.data(s, e, :) = squeeze(dataset(1, b, s, e, x.start : x.end) - dataset(2, b, s, e, x.start : x.end));
            end
        end
        y.data = squeeze(mean(y.data, [1, 3]));

        % select significant electrodes
        chans2plot = [];
        if ~isempty(significant_electrodes)
            for d = 1:length(significant_electrodes)
                for e = 1:length(significant_electrodes(d).times)
                    for f = significant_electrodes(d).times{e}
                        if f >= params.TOI(c, 1) && f <= params.TOI(c, 2)
                            chans2plot(end + 1) = significant_electrodes(d).electrode_n;
                        end
                    end
                end
            end
        end
        
        % plot the topography
        fig = figure(figure_counter);
        topoplot(double(y.data), params.chanlocs, 'style', 'map', 'shading', 'interp', 'whitebk', 'on', ...
                        'maplimits', params.topo_lim, 'electrodes', 'off');
        hold on
        title(sprintf('%s: %s', params.components{c}, params.timepoint{b+1}))
        set(gca, 'FontSize', 14)

        % mark significant electrodes
        if ~isempty(chans2plot)
            for i = 1:length(chans2plot)
                % extract electrode coordinates
                electrode.theta = params.chanlocs(chans2plot(i)).theta - 90;   
                electrode.radius = params.chanlocs(chans2plot(i)).radius;   

                % convert polar to cartesian coordinates
                electrode.x = electrode.radius .* cosd(electrode.theta);  
                electrode.y = -electrode.radius .* sind(electrode.theta);  

                % adjust for topoplot scaling 
                electrode.scale_factor = 0.4 / 0.511;
                electrode.x_scaled = electrode.x * electrode.scale_factor;
                electrode.y_scaled = electrode.y * electrode.scale_factor;

                hold on;
                plot(electrode.x_scaled, electrode.y_scaled, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5);
            end
        end

        % % add colorbar
        % cb = colorbar('Location', 'southoutside');  
        % set(cb, 'Orientation', 'horizontal');  
        % ylabel(cb, 'amplitude (µV)', 'FontSize', 12);

        % save figure and update counter
        saveas(fig, sprintf('%s\\figures\\change_topo_interaction_%s_%s.png', folder.output, params.components{c}, params.timepoint{b+1}))
        figure_counter = figure_counter + 1;
    end
end

% clear and move on
clear a b c d e f i s p header dataset dataset_shuffled toi p_values p_values_correct t_values data_pain data_control perm_cluster_stat perm_clusters ....
        perm perm_p_values perm_cluster_stats combined_samples shuffled_labels cluster significant_samples max_cluster_stats perm_t_values ...
        significant_matrix significant_times chans2plot cluster_indices cluster_t clusters electrode max_perm_cluster_stat neighbors stats ...
        row_counter fig screen_size x y colors idx significant_electrodes ax statement data visual dataset thresholds row_counter_ST
fprintf('section 11 finished.\n')

%% 12) extraction of specific TEP measures
% ----- section input -----
params.prefix = 'icfilt ica ar ffilt sspsir';
params.filename = 'CAPSTEP grand_average';
params.subjects = 20;
% -------------------------
fprintf('section 12: extraction of TEP measures\n')

% update
load(output_file, 'TEP_new_data', 'TEP_new_measures')
cd(folder.output)
clear subject_idx

% load and adjust a header
data2load = dir(sprintf('%s\\%s*%s*.lw6', folder.processed, params.prefix, study));
if ~isempty(data2load)
    load(sprintf('%s\\%s', data2load(1).folder, data2load(1).name), '-mat')
    header.datasize(1) = 1;
    header.events = [];
    params.labels = {header.chanlocs.labels};
    params.x = (0:header.datasize(6)-1)*header.xstep + header.xstart;
else 
    error('ERROR: no suitable headers were wound in the data folder!')
end

% save grand average TEPs for inspection in letswave
fprintf('exporting for letswave: ')
addpath(genpath([folder.toolbox '\letswave 7']));
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        fprintf('. ')

        % export original dataset
        dataset = TEP_new_data.normalized_flipped;
        data(1, :, 1, 1, 1, :) = squeeze(mean(dataset(a, b, :, :, :), 3));
        save(sprintf('%s\\%s %s %s.mat', folder.processed, params.filename, params.condition{a}, params.timepoint{b}), 'data')
        header.name = sprintf('%s %s %s', params.filename, params.condition{a}, params.timepoint{b});
        save(sprintf('%s\\%s %s %s.lw6', folder.processed, params.filename, params.condition{a}, params.timepoint{b}), 'header')

        % export change from baseline
        if b < length(params.timepoint)
            dataset = TEP_new_data.change;
            data(1, :, 1, 1, 1, :) = squeeze(mean(dataset(a, b, :, :, :), 3));
            save(sprintf('%s\\%s change %s %s.mat', folder.processed, params.filename, params.condition{a}, params.timepoint{b+1}), 'data')
            header.name = sprintf('%s change %s %s', params.filename, params.condition{a}, params.timepoint{b+1});
            save(sprintf('%s\\%s change %s %s.lw6', folder.processed, params.filename, params.condition{a}, params.timepoint{b+1}), 'header')
        end
    end
end
fprintf('done.\n')
cd(folder.processed); letswave

% define timepoints, approximative TOIs and EOIs
ok = false; 
measure_counter = 1;
prompt = {'timepoint:' 'time of interest:' 'electrodes of interest:' 'polarity'};
definput = {'t' '' '' 'positive/negative'};
dims = [1 50];
measures = struct([]);
while ~ok
    % encode measure of interest    
    dlgtitle = sprintf('measure %d', measure_counter);    
    input = inputdlg(prompt,dlgtitle,dims,definput);

    % save to output structure
    measures(measure_counter).timepoint = split(input{1}, ',');
    measures(measure_counter).TOI = str2num(input{2});
    measures(measure_counter).EOI = split(input{3}, ',');
    measures(measure_counter).polarity = input{4};

    % ask for continuation
    answer = questdlg('Add another measure?', 'next step',...
        'Add', 'Finish', 'Add');
    switch answer
        case 'Add'
            measure_counter = measure_counter + 1;
        case 'Finish'
            ok = true;
    end
end

% extract information about peak change and precise TOIs
dataset = TEP_new_data.change;
for m = 1:length(measures)
    % identify polarity
    if strcmp(measures(m).polarity, 'positive')
        polarity = 1;
    elseif strcmp(measures(m).polarity, 'negative')
        polarity = -1;
    end

    % calculate grand average change
    data = squeeze(mean(dataset(1, find(strcmp(params.timepoint, measures(m).timepoint)) - 1, :, :, :)));

    % calculate mean change at EOIs
    eois = find(ismember(params.labels, measures(m).EOI));
    measures(m).grand_average = mean(data(eois, :), 1); 

    % extract peak measures
    search_window = find(params.x >= (measures(m).TOI - 0.05) & params.x <= (measures(m).TOI + 0.05));
    [pks, locs] = findpeaks(polarity * measures(m).grand_average(search_window), params.x(search_window), ...
                                           'MinPeakProminence', 0.02, ...   
                                           'WidthReference', 'halfheight');
    if ~isempty(pks)
        % select the most prominent peak
        [~, max_idx] = max(proms); 
        measures(m).peak_latency = locs(max_idx);

        % find where the signal crosses 0
        zero_crossings = find(diff(sign(measures(m).grand_average)) ~= 0); 

        % find all other sufficiently prominent inflection points
        [max_peaks, max_locs] = findpeaks(measures(m).grand_average, params.x, 'MinPeakProminence', 0.02);
        [min_peaks, min_locs] = findpeaks(-measures(m).grand_average, params.x, 'MinPeakProminence', 0.02);
        inflection_points = sort([max_locs, min_locs]);

        % extract onset
        latency_zero = params.x(zero_crossings(find(zero_crossings < find(params.x == locs(max_idx)), 1, 'last')));
        latency_inflection = inflection_points(find(inflection_points < locs(max_idx), 1, 'last'));
        if abs(latency_zero - locs(max_idx)) < abs(latency_inflection - locs(max_idx))
            measures(m).TOI(1) = latency_zero;
        elseif abs(latency_inflection - locs(max_idx)) < abs(latency_zero - locs(max_idx))
            measures(m).TOI(1) = latency_inflection;
        end

        % extract offset
        latency_zero = params.x(zero_crossings(find(zero_crossings > find(params.x == locs(max_idx)), 1, 'first')));
        latency_inflection = inflection_points(find(inflection_points > locs(max_idx), 1, 'first'));
        if abs(latency_zero - locs(max_idx)) < abs(latency_inflection - locs(max_idx))
            measures(m).TOI(2) = latency_zero;
        elseif abs(latency_inflection - locs(max_idx)) < abs(latency_zero - locs(max_idx))
            measures(m).TOI(2) = latency_inflection;
        end

        % plot mean change + TOI
        visual.data = measures(m).grand_average; 
        visual.x = params.x; 
        fig = figure(figure_counter);
        plot_ERP(visual, 'shading', 'off', 'interpolated', measures(m).TOI, ...
            'legend', 'off', 'colours', [0 0 0], 'xlim', [-0.05 0.4])
        title(sprintf('measure %d - %s', m, measures(m).timepoint))

        % save figure and update counter
        saveas(fig, sprintf('%s\\figures\\TEP_measure%d_%s.svg', folder.output, m, measures(m).timepoint))
        figure_counter = figure_counter + 1;
    end
end

% add the very late measure and save
measures(end).TOI = [280, 380];
measures(end).polarity = 'negative';
measures(end).EOI = {header.chanlocs.labels};
TEP_new_measures.TEP_change = measures;  
save(output_file, 'TEP_new_measures', '-append')

% extract single subject mean change of early measures in pain condition
fprintf('extracting early measures: ')
for m = 1:length(measures) - 1 
    % initiate output structure
    component = struct([]);
    fprintf('. ')

    % loop through subjects
    for s = 1:params.subjects
        % encode subject ID 
        component(s).ID = TEP_new(s).ID;

        % select the data at TOI
        b = find(strcmp(params.timepoint, measures(m).timepoint)) - 1;
        eois = find(ismember(params.labels, measures(m).EOI));
        times = find(params.x >= measures(m).TOI(1) & params.x <= measures(m).TOI(2)); 
        data = squeeze(mean(dataset(1, b, s, eois, times), 4))';

        % identify peak latency & amplitude in TOI signal
        if strcmp(measures(m).polarity, 'positive')
            [peak_amplitude, peak_idx] = max(data); 
        elseif strcmp(measures(m).polarity, 'negative')
            [peak_amplitude, peak_idx] = min(data); 
        end
        peak_latency = params.x(find(params.x == measures(m).TOI(1)) + peak_idx);
        % peak_latency = params.x(1526 + peak_idx);

        if ~isempty(peak_latency)            
            % identify timepoints with expected polarity across the full epoch
            data = squeeze(mean(dataset(1, b, s, eois, :), 4))';
            if strcmp(measures(m).polarity, 'positive')
                valid_idx = data > 0; 
            elseif strcmp(measures(m).polarity, 'negative')
                valid_idx = data < 0; 
            end
            onset_idx = find(params.x < peak_latency & ~valid_idx, 1, 'last') + 1;
            offset_idx = find(params.x > peak_latency & ~valid_idx, 1, 'first') - 1;
            component_idx = find(params.x >= params.x(onset_idx) & params.x <= params.x(offset_idx));
            
            % encode to the output structure, if onset before group average 
            % peak latency
            if params.x(onset_idx) <= measures(m).peak_latency
                component(s).peak_latency = peak_latency;
                component(s).peak_amplitude = peak_amplitude;
                component(s).onset = params.x(onset_idx);
                component(s).offset = params.x(offset_idx);
                component(s).duration = params.x(offset_idx) - params.x(onset_idx);
                if component(s).duration > 0
                    component(s).AUC = trapz(params.x(component_idx), data(component_idx));
                else
                    component(s).AUC = NaN;
                end
            else
                % encode NaNs
                component(s).peak_latency = NaN;
                component(s).peak_amplitude = NaN;
                component(s).onset = NaN;
                component(s).offset = NaN;
                component(s).duration = NaN;
                component(s).AUC = NaN;
            end
        else
            % encode NaNs
            component(s).peak_latency = NaN;
            component(s).peak_amplitude = NaN;
            component(s).onset = NaN;
            component(s).offset = NaN;
            component(s).duration = NaN;
            component(s).AUC = NaN;
        end
    end

    % append results to the global output structure
    TEP_new_measures.TEP_change(m).individual_measures = component;

    % calculate group statistics
    subject_idx = ~isnan([component.AUC]);
    TEP_new_measures.TEP_change(m).group_stats.subjects = sum(subject_idx);
    TEP_new_measures.TEP_change(m).group_stats.subject_idx = find(subject_idx);
    TEP_new_measures.TEP_change(m).group_stats.peak_latency.mean = mean([component(subject_idx).peak_latency]);
    TEP_new_measures.TEP_change(m).group_stats.peak_latency.SD = std([component(subject_idx).peak_latency]);
    TEP_new_measures.TEP_change(m).group_stats.peak_latency.SEM = std([component(subject_idx).peak_latency])/sqrt(sum(subject_idx));
    TEP_new_measures.TEP_change(m).group_stats.peak_amplitude.mean = mean([component(subject_idx).peak_amplitude]);
    TEP_new_measures.TEP_change(m).group_stats.peak_amplitude.SD = std([component(subject_idx).peak_amplitude]);
    TEP_new_measures.TEP_change(m).group_stats.peak_amplitude.SEM = std([component(subject_idx).peak_amplitude])/sqrt(sum(subject_idx));
    TEP_new_measures.TEP_change(m).group_stats.onset.mean = mean([component(subject_idx).onset]);
    TEP_new_measures.TEP_change(m).group_stats.onset.SD = std([component(subject_idx).onset]);
    TEP_new_measures.TEP_change(m).group_stats.onset.SEM = std([component(subject_idx).onset])/sqrt(sum(subject_idx));
    TEP_new_measures.TEP_change(m).group_stats.offset.mean = mean([component(subject_idx).offset]);
    TEP_new_measures.TEP_change(m).group_stats.offset.SD = std([component(subject_idx).offset]);
    TEP_new_measures.TEP_change(m).group_stats.offset.SEM = std([component(subject_idx).offset])/sqrt(sum(subject_idx));
    TEP_new_measures.TEP_change(m).group_stats.duration.mean = mean([component(subject_idx).duration]);
    TEP_new_measures.TEP_change(m).group_stats.duration.SD = std([component(subject_idx).duration]);
    TEP_new_measures.TEP_change(m).group_stats.duration.SEM = std([component(subject_idx).duration])/sqrt(sum(subject_idx));
    TEP_new_measures.TEP_change(m).group_stats.AUC.mean = mean([component(subject_idx).AUC]);
    TEP_new_measures.TEP_change(m).group_stats.AUC.SD = std([component(subject_idx).AUC]);
    TEP_new_measures.TEP_change(m).group_stats.AUC.SEM = std([component(subject_idx).AUC])/sqrt(sum(subject_idx));
end
save(output_file, 'TEP_new_measures', '-append')
fprintf('done.\n')

% extract single subject change GFP at the very late measure TOI
dataset = TEP_new_data.change_GFP;
times = find(params.x >= TEP_new_measures.TEP_change(end).TOI(1) & params.x <= TEP_new_measures.TEP_change(end).TOI(2)); 
component = struct([]);
for s = 1:params.subjects
    for a = 1:length(params.condition) 
        for b = 1:length(params.timepoint) - 1
            % extract data
            data = squeeze(dataset(a, b, s, times))';

            % identify peak measures
            [component(a).peak_amplitude(b, s), peak_idx] = max(data);
            component(a).peak_latency(b, s) = params.x(find(params.x == TEP_new_measures.TEP_change(end).TOI(1)) + peak_idx);

            % calculate AUC
            component(a).AUC(b, s) = trapz(params.x(times), data);
        end
    end
end
TEP_new_measures.TEP_change(end).individual_measures = component;

% calculate group stats
row_counter = 1;
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint) - 1    
        TEP_new_measures.TEP_change(end).group_stats(row_counter).condition = params.condition{a};
        TEP_new_measures.TEP_change(end).group_stats(row_counter).timepoint = params.timepoint{b+1};
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_latency.mean = mean([component(a).peak_latency(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_latency.SD = std([component(a).peak_latency(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_latency.SEM = std([component(a).peak_latency(b, :)])/sqrt(20);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_amplitude.mean = mean([component(a).peak_amplitude(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_amplitude.SD = std([component(a).peak_amplitude(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).peak_amplitude.SEM = std([component(a).peak_amplitude(b, :)])/sqrt(20);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).AUC.mean = mean([component(a).AUC(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).AUC.SD = std([component(a).AUC(b, :)]);
        TEP_new_measures.TEP_change(end).group_stats(row_counter).AUC.SEM = std([component(a).AUC(b, :)])/sqrt(20);
        row_counter = row_counter + 1;
    end
end
save(output_file, 'TEP_new_measures', '-append')

% import and adjust MEP results
load(output_file, 'CAPSTEP_MEP')
subject_idx = 1; current_subject = 1;
for c = 1:height(CAPSTEP_MEP)
    if current_subject == str2num(CAPSTEP_MEP.subject{c})
        CAPSTEP_MEP.ID{c} = TEP_new(subject_idx).ID;
    elseif current_subject < str2num(CAPSTEP_MEP.subject{c})
        current_subject = str2num(CAPSTEP_MEP.subject{c});
        subject_idx = subject_idx + 1;
        CAPSTEP_MEP.ID{c} = TEP_new(subject_idx).ID;
    end
end

% export for R together with MEPs 
TEP_new_export = table;
row_counter = height(TEP_new_export) + 1;
for m = 1:length(TEP_new_measures.TEP_change) 
    if m < length(TEP_new_measures.TEP_change) 
        for s = 1:params.subjects
            if ismember(s, TEP_new_measures.TEP_change(m).group_stats.subject_idx)
                % encode subject info
                TEP_new_export.measure(row_counter) = m;
                TEP_new_export.subject{row_counter} = TEP_new(s).ID;
                TEP_new_export.condition{row_counter} = 'pain';
                TEP_new_export.timepoint{row_counter} = TEP_new_measures.TEP_change(m).timepoint;

                % encode TEPs
                TEP_new_export.peak_latency(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).peak_latency;
                TEP_new_export.peak_amplitude(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).peak_amplitude;
                TEP_new_export.onset(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).onset;
                TEP_new_export.offset(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).offset;
                TEP_new_export.duration(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).duration;
                TEP_new_export.AUC(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(s).AUC;

                % encode MEPs
                TEP_new_export.MEP_change(row_counter) = CAPSTEP_MEP.amplitude_norm(strcmp(CAPSTEP_MEP.ID, TEP_new(s).ID) ...
                    & strcmp(CAPSTEP_MEP.session, 'caps') ...
                    & strcmp(CAPSTEP_MEP.timepoint, TEP_new_measures.TEP_change(m).timepoint)) / 100;
                
                % update counter
                row_counter = row_counter + 1;
            end
        end
    else
        for s = 1:params.subjects
            for a = 1:length(params.condition)
                % identidy session label
                if a == 1
                    session = 'caps';
                elseif a == 2
                    session = 'ctrl';
                end

                for b = 1:length(params.timepoint) - 1 
                    % encode subject info
                    TEP_new_export.measure(row_counter) = m;
                    TEP_new_export.subject{row_counter} = TEP_new(s).ID;
                    TEP_new_export.condition{row_counter} = params.condition{a};
                    TEP_new_export.timepoint{row_counter} = params.timepoint{b+1};

                    % encode TEPs
                    TEP_new_export.peak_latency(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(a).peak_latency(b, s);
                    TEP_new_export.peak_amplitude(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(a).peak_amplitude(b, s);
                    TEP_new_export.AUC(row_counter) = TEP_new_measures.TEP_change(m).individual_measures(a).AUC(b, s);
                    
                    % encode MEPs
                    TEP_new_export.MEP_change(row_counter) = CAPSTEP_MEP.amplitude_norm(strcmp(CAPSTEP_MEP.ID, TEP_new(s).ID) ...
                        & strcmp(CAPSTEP_MEP.session, session) ...
                        & strcmp(CAPSTEP_MEP.timepoint, params.timepoint{b+1})) / 100;
                    
                    % update counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
writetable(TEP_new_export, sprintf('%s\\TEP_new_export.csv', folder.output))
save(output_file, 'TEP_new_export', '-append')

% clear and move on
clear a b c m s dataset data2load header data ok prompt dlgtitle dims definput input measures measure_counter ...
    answer current_subject CAPSTEP_MEP times search_window max_idx fig latency_zero latency_inflection session...
    polarity eois pks locs zero_crossings max_peaks max_locs min_peaks min_locs inflection_points ...
    peak_amplitude peak_idx peak_latency valid_idx onset_idx offset_idx component_idx component subject_idx AUC visual row_counter
fprintf('section 12 finished.\n') 

%% code scraps
% adjust history for ICA
for a = 2:length(params.condition) 
    fprintf('%s session: ', params.condition{a})
    for b = 1:length(params.timepoint)
        fprintf('. ')
        load(sprintf('ica ar ffilt sspsir CAPSTEP S06 %s %s .lw6', params.condition{a}, params.timepoint{b}), '-mat')
        header.history(10).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
        header.history(10).configuration.parameters = header.history(10).option;  
        header.history(10).configuration.parameters.ICA_um = header.history(10).configuration.parameters.unmix_matrix; 
        header.history(10).configuration.parameters.ICA_mm = header.history(10).configuration.parameters.mix_matrix; 
        header.history(10).configuration.parameters = rmfield(header.history(10).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
        save(sprintf('ica ar ffilt sspsir CAPSTEP S06 %s %s .lw6', params.condition{a}, params.timepoint{b}), 'header');
    end
    fprintf('done.\n')
end

% adjust eventcodes
for a = 1:length(params.condition) 
    for b = 1:length(dataset(a).raw)
        for e = 1:length(dataset(a).raw(b).header.events)
            dataset(a).raw(b).header.events(e).code = 'TMS'; 
        end
    end
end

% remove bad epochs in EEGLAB
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        pop_eegplot(EEG, 1, 1, 1);
        waitfor(gcf); 
        EEG = eeg_checkset(EEG);
        
        % encode bad epochs to info structure
        answer = questdlg('Have you discarded any epochs?', 'Bad epochs',...
            'YES', 'NO', 'NO');
        if a == 1 && b == 1
            TEP_new(subject_idx).processing(13).process = sprintf('bad trials discarded');
            TEP_new(subject_idx).processing(13).params.GUI = 'EEGLAB';
            TEP_new(subject_idx).processing(13).suffix = params.suffix{4};
            TEP_new(subject_idx).processing(13).date = sprintf('%s', date);
        end
        switch answer
            case 'YES'
                for c = 1:length(ALLCOM)
                    if contains(ALLCOM{c}, 'pop_rejepoch')
                        match = regexp(ALLCOM{c}, '\[(.*?)\]', 'match');
                        if ~isempty(match)
                            discarded = str2num(match{1}(2:end-1)); 
                        else
                            discarded = []; 
                        end
                        break
                    end
                end
            case 'NO'
                discarded = []; 
        end
        
        TEP_new(subject_idx).processing(13).params.discarded{a, b} = discarded;
    end
end

% encode bad trials
check = true;
while check
    % ask for bad trials
    if exist('input_old') ~= 1
        for a = 1:length(params.condition)
            for b = 1:length(params.timepoint)
                prompt{(a-1)*length(params.timepoint) + b} = sprintf('%s - %s', params.condition{a}, params.timepoint{b});
                definput{(a-1)*length(params.timepoint) + b} = '';
            end
        end
        dlgtitle = 'bad trials';
        dims = [1 60];
        input = inputdlg(prompt,dlgtitle,dims,definput);
    else
        for a = 1:length(params.condition)
            for b = 1:length(params.timepoint)
                prompt{(a-1)*length(params.timepoint) + b} = sprintf('%s - %s', params.condition{a}, params.timepoint{b});
                definput{(a-1)*length(params.timepoint) + b} = input_old{(a-1)*length(params.timepoint) + b};
            end
        end
        dlgtitle = 'bad trials';
        dims = [1 60];
        input = inputdlg(prompt,dlgtitle,dims,definput);
    end

    % verify if the numbers match
    for a = 1:length(params.condition)
        for b = 1:length(params.timepoint)
            discarded = str2num(input{(a-1)*length(params.timepoint) + b});
            if size(dataset(a).sspsir(b).data, 1) - length(discarded) == size(dataset(a).checked(b).data, 1)
                match(a, b) = 0; 
            else
                match(a, b) = 1; 
            end
        end
    end
    
    % encode if match
    if sum(match, 'all') == 0
        % encode 
        fprintf('encoding bad trials... ')
        TEP_new(subject_idx).processing(13).process = sprintf('bad trials discarded');
        TEP_new(subject_idx).processing(13).params.GUI = 'letswave';
        TEP_new(subject_idx).processing(13).suffix = params.suffix{1};
        TEP_new(subject_idx).processing(13).date = sprintf('%s', date);
        for a = 1:length(params.condition)
            for b = 1:length(params.timepoint)
                TEP_new(subject_idx).processing(13).params.discarded{a, b} = str2num(input{(a-1)*length(params.timepoint) + b});
                TEP_new(subject_idx).processing(13).params.kept(a, b) = size(dataset(a).checked(b).data, 1);
            end
        end

        % exit the loop
        check = false;
    else
        % provide information
        fprintf('ATTENTION - in following datasets the number of retained epochs does not match the number of discarded expochs:\n')
        for a = 1:length(params.condition)
            for b = 1:length(params.timepoint)
                if match(a, b) == 1
                    fprintf('%s - %s\n', params.condition{a}, params.timepoint{b})
                end
            end
        end
        fprintf('please correct your input\n')

        % store previous responses
        input_old = input;
    end
end

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
    if length(datas) == length(headers) 
        for d = 1:length(datas)
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
        EEG.event(e).latency = (e-1)*EEG.pnts + EEG.xmin*(-1)*EEG.srate;
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
function lwdata = export_lw(EEG, header, name)
% =========================================================================
% exports data from EEGLAB to letswave format
% ========================================================================= 
% transform data
data = [];
for t = 1:size(EEG.data, 3)
    for e = 1:size(EEG.data, 1)
        for i = 1:size(EEG.data, 2)
            data(t, e, 1, 1, 1, i) = EEG.data(e, i, t);
        end
    end
end
lwdata.data = data;    
    
% modify header
lwdata.header = header; 
lwdata.header.name = name;
lwdata.header.datasize = size(lwdata.data);
lwdata.header.chanlocs = lwdata.header.chanlocs(1:size(lwdata.data, 2));
lwdata.header.events = lwdata.header.events(1:size(lwdata.data, 1));
end
function plot_ERP(input, varargin)
% =========================================================================
% plots an event-related potential
% input = structure with fields:    
%           data --> condition/electrode * sample
%           x --> vector with time samples
%           CI_upper --> condition/electrode * sample
%           CI_lower --> condition/electrode * sample
% varargins = name-value pairs: 
%           xlim --> 2-element vector (min, max)     
%           ylim --> 2-element vector (min, max) 
%           colours --> n*3 matrix of RGB values
%           shading --> 'on'(default)/'off'
%           alpha --> a float (default 0.2)           
%           legend --> 'on'(default)/'off'
%           labels --> cell array with labels for the legend  
%           legend_loc --> legend location (default 'southeast')
%           eoi --> label of a channel to be highlighted
%           reverse --> 'on'/'off'(default) - flips y axis
%           interpolated --> time window that was interpolated
% =========================================================================  
% set defaults
x_limits = [0,0];
y_limits = [0,0];
col = prism(size(input.data, 1));
shading = true;
alpha = 0.2;
plot_legend = true;
for c = 1:size(input.data, 1)
    labels{c} = sprintf('condition %d', c);
end
legend_loc = 'southeast';
highlight = false;
reverse = false;
interpolate = false;

% check for varargins
if ~isempty(varargin)
    % x limits
    a = find(strcmpi(varargin, 'xlim'));
    if ~isempty(a)
        x_limits = varargin{a + 1};
    end

    % y limits
    b = find(strcmpi(varargin, 'ylim'));
    if ~isempty(b)
        y_limits = varargin{b + 1};
    end

    % colours
    c = find(strcmpi(varargin, 'colours'));
    if ~isempty(c)
        col = varargin{c + 1};
    end

    % shading - default on
    d = find(strcmpi(varargin, 'shading'));
    if ~isempty(d) && strcmp(varargin{d + 1}, 'off')
        shading = false;
    end

    % alpha
    e = find(strcmpi(varargin, 'alpha'));
    if ~isempty(e)
        alpha = varargin{e + 1};
    end

    % legend - default on
    f = find(strcmpi(varargin, 'legend'));
    if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
        plot_legend = false;
    end    

    % labels
    g = find(strcmpi(varargin, 'labels'));
    if ~isempty(g)
        labels = varargin{g + 1};
    end

    % legend location
    h = find(strcmpi(varargin, 'legend_loc'));
    if ~isempty(h) 
        legend_loc = varargin{h + 1};
    end  

    % highlighted channel - default off
    i = find(strcmpi(varargin, 'eoi'));
    if ~isempty(i)
        eoi = varargin{i + 1};
        eoi_n = find(contains(input.chanlabels, eoi));
        highlight = true;
    end 

    % interpolated interval - default off
    j = find(strcmpi(varargin, 'interpolated'));
    if ~isempty(j)
        interpolate_toi = varargin{j + 1};
        interpolate = true;
    end 

    % reverse y axis - default off
    r = find(strcmpi(varargin, 'reverse'));
    if ~isempty(r) && strcmp(varargin{r + 1}, 'on')
        reverse = true;
    end
end

% loop through datasets to plot
for t = 1:size(input.data, 1) 
    P(t) = plot(input.x, input.data(t, :), 'Color', col(t, :), 'LineWidth', 2);
    hold on
    if shading
        F(t) = fill([input.x fliplr(input.x)],[input.CI_upper(t, :) fliplr(input.CI_lower(t, :))], ...
            col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
end

% check y limits
if y_limits(1) == 0 && y_limits(2) == 0
    y_limits = ylim;
end

% highlight channel if required
if highlight
    P(end + 1) = plot(input.x, input.data(eoi_n, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 4);
    hold on
end

% shade interpolated window if required
if interpolate
    interpolate_x = [interpolate_toi(1), interpolate_toi(2), interpolate_toi(2), interpolate_toi(1)];
    interpolate_y = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
    fill(interpolate_x, interpolate_y, [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')

% plot legend if required
if plot_legend 
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
else
    legend('off')
end

% axes
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.TickDir = 'out'; 
ax.XColor = [0.5020    0.5020    0.5020]; 
ax.YColor = [0.5020    0.5020    0.5020]; 

% set x limits 
if x_limits(1) == 0 && x_limits(2) == 0
    xlim([input.x(1), input.x(end)]) 
else
    xlim(x_limits)
end

% referse y axis if required
if reverse
    set(gca, 'YDir', 'reverse');
end

% other parameters
xlabel('time (s)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
set(gca, 'Layer', 'Top')
end
function clusters = find_contiguous_clusters(significant_samples)
% =========================================================================
% identify clusters of consecutive samples
% significant_samples = a vector of indexes of samples with p < alpha
% =========================================================================  
if isempty(significant_samples)
    return; 
end    
clusters = {};
current_cluster = [significant_samples(1)];
for i = 2:length(significant_samples)
    if significant_samples(i) == significant_samples(i-1) + 1
        current_cluster = [current_cluster, significant_samples(i)];
    else
        clusters{end+1} = current_cluster;
        current_cluster = [significant_samples(i)];
    end
end
clusters{end+1} = current_cluster; 
end
function clusters = find_spatiotemporal_clusters(significant_matrix, neighbors)
% =========================================================================
% identify connected clusters in time and across electrodes
% input:    significant_matrix = binary matrix (electrodes x timepoints) where 1 = significant  
%           neighbors = cell array with neighbor indices for each electrode
% =========================================================================  
clusters = {};
[elecs, times] = size(significant_matrix);
checked = zeros(size(significant_matrix)); 
for e = 1:elecs
    for t = 1:times
        if significant_matrix(e, t) == 1 && checked(e, t) == 0
            % start a new cluster
            cluster = [];
            [cluster, checked] = grow_cluster(significant_matrix, e, t, cluster, checked, neighbors);
            clusters{end+1} = cluster; 
        end
    end
end
end
function [cluster, checked] = grow_cluster(significant_matrix, e, t, cluster, checked, neighbors)
% =========================================================================
% recursive function to grow clusters
% input:    significant_matrix = binary matrix (electrodes x timepoints) where 1 = significant  
%           e = current electrode
%           t = current timepoint
%           cluster = active cluster
%           checked = bin checked = 1/not checked = 0
%           neighbors = cell array with neighbor indices for each electrode
% =========================================================================
[elecs, times] = size(significant_matrix);   
if e < 1 || e > elecs || t < 1 || t > times || significant_matrix(e, t) == 0 || checked(e, t) == 1
    return;
end
    
% mark as checked and add to cluster
checked(e, t) = 1;
cluster = [cluster; e, t];

% check time neighbors
time_neighbors = [t-1, t+1];
for tn = 1:length(time_neighbors)
    if time_neighbors(tn) >= 1 && time_neighbors(tn) <= times
        [cluster, checked] = grow_cluster(significant_matrix, e, time_neighbors(tn), cluster, checked, neighbors);
    end
end

% check spatial eighbors
for n = neighbors{e} 
    [cluster, checked] = grow_cluster(significant_matrix, n, t, cluster, checked, neighbors);
end

end
function neighbors = find_electrode_neighbors(chanlocs, max_distance)
% =========================================================================
% identifies electrode neighbors based on spatial proximity.
% max_distance = how close electrodes must be to be considered neighbors.
% =========================================================================  
num_electrodes = length(chanlocs);
neighbors = cell(num_electrodes, 1);
for i = 1:num_electrodes
    for j = 1:num_electrodes
        % compute Euclidean distance between electrodes
        dist = sqrt((chanlocs(i).X - chanlocs(j).X)^2 + ...
                    (chanlocs(i).Y - chanlocs(j).Y)^2 + ...
                    (chanlocs(i).Z - chanlocs(j).Z)^2);

        % if within threshold and not itself, store as neighbor
        if dist > 0 && dist < max_distance
            neighbors{i} = [neighbors{i}, j];
        end
    end
end
end
function plot_distribution(visual)
% =========================================================================
% plots a null distribution of values obtained via permutations, marks the
% significance threshold
% input:    visual = structure with fields: 
%                       - values
%                       - threshold
%                       - color
% =========================================================================  
% plot the histogram of real values
histogram(visual.values, 'Normalization', 'pdf', 'FaceColor', visual.color, 'FaceAlpha', 0.4, 'EdgeColor', [1 1 1]); 

% fit a smooth curve to the distribution
x_range = linspace(min(visual.values), max(visual.values), 100);
pdf_fit = ksdensity(visual.values, x_range); % kernel density estimate
plot(x_range, pdf_fit, 'Color', visual.color, 'LineWidth', 2);

% mark the threshold
x_step = (max(visual.values) - min(visual.values))/100; 
xline(visual.threshold, 'k--', 'LineWidth', 2);
text(visual.threshold + x_step*8, max(pdf_fit) * 0.9, sprintf('%.3f', visual.threshold), ...
'Color', 'k', 'FontSize', 12, 'Rotation', 90, 'FontWeight', 'bold');

% labels + title
xlabel('cluster size distribution');
ylabel('probability density');
title(visual.title);

% appearance
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.TickDir = 'out'; 
ax.XColor = [0.5020    0.5020    0.5020]; 
ax.YColor = [0.5020    0.5020    0.5020]; 
set(gca, 'FontSize', 12)
set(gca, 'Layer', 'Top')
end