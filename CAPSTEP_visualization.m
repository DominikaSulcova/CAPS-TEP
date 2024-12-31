%% CAPS-TEP visualizazion

%% parameters
% load the output files
load('CAPSTEP_output.mat', 'CAPSTEP_TEP_data', 'CAPSTEP_TEP_default', 'CAPSTEP_MEP_data')

% directory
folder.data = uigetdir(pwd, 'Choose the input folder');
folder.figures = sprintf('%s\\figures', folder.data);

% dataset
params.participant = 1:20;
params.condition = {'caps' 'ctrl'};
params.timepoint = {'baseline' 't1' 't2' 't3' 't4' 't5' 't6'};
params.peaks = CAPSTEP_TEP_default.peaks;
params.electrodes = CAPSTEP_TEP_default.electrodes;

% visualization
visual.eoi = 'Cz';
visual.time_window = [-0.05, 0.3];
visual.xstep = 0.0005;
visual.x = [visual.time_window(1) : visual.xstep : visual.time_window(2)];
visual.z = 1.96;
visual.alpha = 0.15;
visual.figure_counter = 1;
visual.colours = [0.8902    0.2118    0.2118; 0    0.4471    0.7412];

%% all-average TEP
% select average data
visual.data = squeeze(mean(CAPSTEP_TEP_data, [1:3]));

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
plot_ERP(visual.data, [], [], visual.x, 'colours', visual.colour_butterfly, 'legend', 'off', 'eoi', visual.eoi_n)

% save figure and update counter
saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
visual.figure_counter = visual.figure_counter + 1; 
clear c fig figname 

%% all-average GFP
% select average data
visual.data = squeeze(mean(CAPSTEP_TEP_data, [1:3]));
visual.data = squeeze(std(visual.data, 1));

% launch the figure
figname = 'CAPSTEP_TEP_all_GFP';
fig = figure(visual.figure_counter);
set(fig, 'Name', figname, 'Position', [100, 100, 800, 550])

% plot butterfly plot
plot_ERP(visual.data, [], [], visual.x, 'colours', [0 0 0], 'legend', 'off')

% save figure and update counter
saveas(fig, sprintf('%s\\%s.svg', folder.figures, figname))
visual.figure_counter = visual.figure_counter + 1; 
clear c fig figname 

%% plot timecourse
% loop through each timepoint
for t = 1:length(time)
    % select data
    coi = find(strcmp(electrodes, eoi));
    t_value = tinv(0.975, size(data, 2) - 1); 
    data = squeeze(CAPSTEP_TEP_data(:, t, :, coi, :));
    visual_data = squeeze(mean(data, 2));
    visual_sem(1, :) = visual_data(1, :) + std(squeeze(data(1, :, :)), 0, 1) / sqrt(size(data, 2)); 
    visual_sem(2, :) = visual_data(2, :) + std(squeeze(data(2, :, :)), 0, 1) / sqrt(size(data, 2)); 
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
function plot_ERP(visual_data, visual_CI_upper, visual_CI_lower, x, varargin)
% =========================================================================
% plots an event-related potential
% =========================================================================  
% set defaults
col = prism(size(visual_data, 1));
alpha = 0.2;
y_limits = [0,0];
for c = 1:size(visual_data, 1)
    labels{c} = sprintf('condition %d', c);
end
shading = true;
plot_legend = true;
legend_loc = 'southeast';
highlight = false;

% check for varargins
if ~isempty(varargin)
    % y limits
    a = find(strcmpi(varargin, 'ylim'));
    if ~isempty(a)
        y_limits = varargin{a + 1};
    end

    % labels
    b = find(strcmpi(varargin, 'labels'));
    if ~isempty(b)
        labels = varargin{b + 1};
    end

    % colours
    d = find(strcmpi(varargin, 'colours'));
    if ~isempty(d)
        col = varargin{d + 1};
    end

    % alpha
    e = find(strcmpi(varargin, 'alpha'));
    if ~isempty(e)
        alpha = varargin{e + 1};
    end

    % shading - default on
    f = find(strcmpi(varargin, 'shading'));
    if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
        shading = false;
    end

    % legend - default on
    f = find(strcmpi(varargin, 'legend'));
    if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
        plot_legend = false;
    end       

    % legend location
    g = find(strcmpi(varargin, 'legend_loc'));
    if ~isempty(g) 
        legend_loc = varargin{g + 1};
    end  
    % highlighted channel
    h = find(strcmpi(varargin, 'eoi'));
    if ~isempty(h)
        eoi = varargin{h + 1};
        highlight = true;
    end 
end

% in case there are no CIs
if isempty(visual_CI_upper) || isempty(visual_CI_lower)
    shading = false;
end

% loop through datasets to plot
for t = 1:size(visual_data, 1) 
    P(t) = plot(x, visual_data(t, :), 'Color', col(t, :), 'LineWidth', 2);
    hold on
    if shading
        F(t) = fill([x fliplr(x)],[visual_CI_upper(t, :) fliplr(visual_CI_lower(t, :))], ...
        col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
end

% check y limits
if y_limits(1) == 0 && y_limits(2) == 0
    y_limits = ylim;
end

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--')

% highlight channel if required
if highlight
    P(end + 1) = plot(x, visual_data(eoi, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 3);
end

% legend
if plot_legend
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
end

% axes
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
% ax.TickDir = 'out'; 
ax.XColor = [0.5020    0.5020    0.5020]; 
ax.YColor = [0.5020    0.5020    0.5020]; 

% other parameters
xlabel('time (ms)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
xlim([x(1), x(end)])  
set(gca, 'Layer', 'Top')
% set(gca, 'YDir', 'reverse');
end