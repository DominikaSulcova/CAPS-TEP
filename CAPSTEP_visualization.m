%% CAPS-TEP visualizazion

%% parameters
% load the output files
load('CAPSTEP_output.mat', 'CAPSTEP_TEP_data', 'CAPSTEP_TEP_default', 'CAPSTEP_MEP_data')

% directory
folder.data = uigetdir(pwd, 'Choose the input folder');
folder.figures = sprintf('%s\\figures', folder.data);

% dataset
params.participant = 1:20;
params.condition = {'capsaicin' 'vehicle'};
params.timepoint = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
params.peaks = CAPSTEP_TEP_default.peaks;
params.electrodes = CAPSTEP_TEP_default.electrodes;

% visualization
visual.eoi = 'Cz';
visual.toi_TEP = [-0.05, 0.3];
visual.xstep_TEP = 0.0005;
visual.toi_MEP = [-0.03, 0.1];
visual.xstep_MEP = 0.001;
visual.z = 1.96;
visual.alpha = 0.15;
visual.figure_counter = 1;
visual.colours = [0.8902    0.2118    0.2118; 0    0.4471    0.7412];

%% all-average TEP
% select average data and x
visual.data = squeeze(mean(CAPSTEP_TEP_data, [1:3]));
visual.x = visual.toi_TEP(1) : visual.xstep_TEP : visual.toi_TEP(2);

% choose colours
for c = 1:size(visual.data, 1)
    visual.colour_butterfly(c, :) = [0.3020    0.7451    0.9333];
end

% identify eoi number
visual.eoi_n = find(strcmp(params.electrodes, visual.eoi));

% launch the figure
figname = 'CAPSTEP_TEP_all_avg';
fig = figure(visual.figure_counter);
set(fig, 'Name', figname, 'Position', [100, 100, 800, 550])

% plot butterfly plot
plot_ERP(visual, 'colours', visual.colour_butterfly, 'shading', 'off', 'legend', 'off', 'eoi', visual.eoi_n)

% save figure and update counter
saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
visual.figure_counter = visual.figure_counter + 1; 
clear c fig figname 

%% all-average GFP
% select average data and x
visual.data = squeeze(mean(CAPSTEP_TEP_data, [1:3]));
visual.data = squeeze(std(visual.data, 1));
visual.x = visual.toi_TEP(1) : visual.xstep_TEP : visual.toi_TEP(2);

% launch the figure
figname = 'CAPSTEP_TEP_all_GFP';
fig = figure(visual.figure_counter);
set(fig, 'Name', figname, 'Position', [100, 100, 800, 550])

% plot butterfly plot
plot_ERP(visual, 'colours', [0 0 0], 'shading', 'off', 'legend', 'off')

% save figure and update counter
saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
visual.figure_counter = visual.figure_counter + 1; 
clear c fig figname 

%% timepoint-average MEP
% select data and flip subject n.1
data.original = CAPSTEP_MEP_data;
data.original(2, :, 1, :) = data.original(2, :, 1, :)*-1;

% calculate t value
t_value = tinv(0.975, size(data.original, 3) - 1); 

% calculate mean and CI
for a = 1:size(data.original, 1)
    for b = 1:size(data.original, 2)
        data.mean(a, b, :) = squeeze(mean(data.original(a, b, :, :), 3))';
        data.sem(a, b, :) = squeeze(std(data.original(a, b, :, :), 0, 3)) / sqrt(size(data.original, 3))'; 
        data.CI_upper(a, b, :) = data.mean(a, b, :) + t_value * data.sem(a, b, :); 
        data.CI_lower(a, b, :) = data.mean(a, b, :) - t_value * data.sem(a, b, :); 
    end
end

% select x
visual.x = visual.toi_MEP(1) : visual.xstep_MEP : visual.toi_MEP(2) - visual.xstep_MEP;

% launch the figure 
figname = 'CAPSTEP_TEP_MEP';
fig = figure(visual.figure_counter);
set(fig, 'Name', figname, 'Position', [100, 100, 1000, 550])

% plot
x = -0.2 : visual.xstep_MEP : 0.107;
x_idx = x >= visual.toi_MEP(1) & x < visual.toi_MEP(2);
y_limits = [];
for a = 1:size(data.original, 1)
    for b = 1:size(data.original, 2)
        % subset data
        visual.data = squeeze(data.mean(a, b, x_idx))';
        visual.CI_upper = squeeze(data.CI_upper(a, b, x_idx))';
        visual.CI_lower = squeeze(data.CI_lower(a, b, x_idx))';

        % plot timepoint
        subplot(size(data.original, 1), size(data.original, 2), (a-1)*size(data.original, 2) + b)
        plot_ERP(visual, 'colours', visual.colours(a, :), 'legend', 'off')
        title(sprintf('%s - %s', params.condition{a}, params.timepoint{b}))

        % encode y limits
        y_limits(end+1, :) = get(gca, 'YLim');
    end
end

% set scaling across plots
for a = 1:size(data.original, 1)
    for b = 1:size(data.original, 2)
        subplot(size(data.original, 1), size(data.original, 2), (a-1)*size(data.original, 2) + b)       
        ylim([min(y_limits(:, 1)), max(y_limits(:, 2))])            
    end
end

% save figure and update counter
saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
visual.figure_counter = visual.figure_counter + 1; 
clear a b data t_value figname fig x x_idx y_limits 

%% plot timecourse
% loop through each timepoint
for t = 1:length(time)
    % select data
    coi = find(strcmp(electrodes, eoi));
    t_value = tinv(0.975, size(data, 2) - 1); 
    data = squeeze(CAPSTEP_TEP_data(:, t, :, coi, :));
    visual_data = squeeze(mean(data, 2));
    visual_sem(1, :) = std(squeeze(data(1, :, :)), 0, 1) / sqrt(size(data, 2)); 
    visual_sem(2, :) = std(squeeze(data(2, :, :)), 0, 1) / sqrt(size(data, 2)); 
    visual_CI_upper(1, :) = visual_data(1, :) + t_value * visual_sem(1, :); 
    visual_CI_lower(1, :) = visual_data(1, :) - t_value * visual_sem(1, :); 
    visual_CI_upper(2, :) = visual_data(2, :) + t_value * visual_sem(2, :); 
    visual_CI_lower(2, :) = visual_data(2, :) - t_value * visual_sem(2, :); 

    % launch the figure
    figname = sprintf('CAPSTEP_TEP_timecourse_%s', time{t});
    fig = figure(figure_counter);
    set(fig, 'Name', figname, 'Position', [100, 100, 800, 550])

    % plot the average data
    plot_ERP(visual_data, visual_CI_upper, visual_CI_lower, x, 'colours', colours, 'labels', condition, 'shading', 'off')

    % save figure adn update counter
    saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
    figure_counter = figure_counter + 1; 
end
%% plot peak topographies
CAPSTEP_TEP_peaks.latency 

%% functions
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
%           plot_legend --> 'on'(default)/'off'
%           labels --> cell array with labels for the legend  
%           legend_loc --> legend location (default 'southeast')
%           eoi --> label of a channel to be highlighted
%           reverse --> 'on'/'off'(default) - flips y axis
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
        highlight = true;
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

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')

% highlight channel if required
if highlight
    P(end + 1) = plot(input.x, input.data(eoi, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 3);
end

% plot legend if required
if plot_legend 
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
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
xlabel('time (ms)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
set(gca, 'Layer', 'Top')
end