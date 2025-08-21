
clc
clearvars
close all



%% ------------------ User/session parameters --------------------------------------------
% Basic metadata
animal_list =  {'hans'};        % {'hans', 'chip'}
area_list = {'LIP'};            % {'LIP', 'VIP', 'STS'}
% animal_list =  {'chip'};        % {'hans', 'chip'}
% area_list = {'LIP'};            % {'LIP', 'VIP', 'STS'}

% Analysis tags
settings.force_extraction = 0;
settings.verbosity = 1;         % 0=none, 1=some, 2=lots

% Sampling rates [Hz]
settings.sr_eye = 500;          % eye sampling
settings.sr_neu = 1000;         % neural sampling (samples per second)

% Eye-tracking settings
settings.eyetracking.eye = 3;
settings.eyetracking.speed_threshold_high_degsec = 200;
settings.eyetracking.speed_threshold_low_degsec = 50;
settings.eyetracking.filter_width = 2;
settings.eyetracking.saccade_window_min = -100;
settings.eyetracking.saccade_window_max = +200;
settings.eyetracking.position_window_ms = 50;

% Stimulus channels priority (try first, then fallback)
settings.stim_channels = [7, 4];
settings.min_triggers_delay_ms = 500;

% Spikes parameters
 settings.spikes_filter_width = 10;

% Conditions of interest and their spatial positions (deg)
settings.condition_list = [1:6, 13];
settings.stim_positions = [-25, -15, -5, +5, +15, +25];

% Binning and analysis window
% settings.bin_size_ms = 10;           % bin width [ms]
% settings.bin_size_samples = round(settings.bin_size_ms * settings.sr_neu / 1000); % [samples per bin]

% Analysis windows
settings.time_window    = [-500, 500];     % analysis window relative to alignment [ms]
settings.n_bins         = settings.time_window(2)-settings.time_window(1);
settings.time_vec       = settings.time_window(1):settings.time_window(2)-1; % bin edges [ms]
settings.baseline_window = [-100,-50];
settings.pre_saccadic_window = [-75,-25];
settings.saccadic_window = [-25,+25];
settings.post_saccadic_window = [+25,+75];


% Stimulus alignment delay (applied only to 'stimulus' alignment)
settings.delay_ms      = 200;                                   % [ms]
settings.delay_samples = round(settings.delay_ms * settings.sr_neu / 1000);       % [samples]
settings.perisacc_thresh_ms = 200;                              % [ms]

% Firing rate threshold for active neuron mask
settings.min_FR = 3;  % [Hz]

% Plotting options
settings.plot.n_rows = 4;
settings.plot.n_cols = 4;
settings.fr_max = 80;

% show_summary = true; % se false → nessuna stampa
% save_summary_to_file = true;
% summary_filename = 'pipeline_summary.txt';


% plot_data = 0;                  % plot figures
%     response_window = 200;          % in ms
%     saccade_baseline_window = [-300,+200];
%     saccade_response_window = [-200,+200];
%     visual_window = [-200,+200];
%     baseline_window = [-200,-100];
%     stim_position = [-25:10:+25];   % in deg
%     sacc_dist_threshold = 4;
%     fr_threshold = 4.0;
%     fr_max = 80;
%     area_colors = [sscanf('1b9e77', '%2x%2x%2x'),...
%         sscanf('d95f02', '%2x%2x%2x'),...
%         sscanf('7570b3', '%2x%2x%2x'),...
%         sscanf('e7298a', '%2x%2x%2x')]'/255;
%     location_colors = [sscanf('440154', '%2x%2x%2x'),...
%         sscanf('414487', '%2x%2x%2x'),...
%         sscanf('2a788e', '%2x%2x%2x'),...
%         sscanf('22a884', '%2x%2x%2x'),...
%         sscanf('7ad151', '%2x%2x%2x'),...
%         sscanf('fde725', '%2x%2x%2x')];

% Brewer colormap
area_colors = [sscanf('1b9e77', '%2x%2x%2x'),...
    sscanf('d95f02', '%2x%2x%2x'),...
    sscanf('7570b3', '%2x%2x%2x'),...
    sscanf('e7298a', '%2x%2x%2x')]'/255;

% Viridis colormap
location_colors = [...
    sscanf('440154', '%2x%2x%2x'),...
    sscanf('414487', '%2x%2x%2x'),...
    sscanf('2a788e', '%2x%2x%2x'),...
    sscanf('22a884', '%2x%2x%2x'),...
    sscanf('7ad151', '%2x%2x%2x'),...
    sscanf('fde725', '%2x%2x%2x')]/255;





%% STEP 1: data extraction
fprintf('\n________________\nSTEP 1\nLoad datasets\n');
dataset_num = 0;
dataset_name = cell(1, length(animal_list)*length(area_list));
datasets = cell(1, length(animal_list)*length(area_list));
collapsed = cell(1, length(animal_list)*length(area_list));

for animal_num=1:length(animal_list)

    for area_num=1:length(area_list)
        dataset_num = dataset_num+1;
        dataset_name{dataset_num} = sprintf('dataset_%s_%s.mat', animal_list{animal_num}, area_list{area_num});
        fprintf('%i/%i\n', dataset_num, length(animal_list)*length(area_list));

        if ~exist(['../', dataset_name{dataset_num}], 'file') || settings.force_extraction
            fprintf('Extracting data...\n');
            dataset = extract_dataset(animal_list{animal_num}, area_list{area_num}, settings);
            save(['../', dataset_name{dataset_num}], 'dataset', '-v7.3');
            datasets{dataset_num} = dataset;
            clear dataset
        else
            load(['../', dataset_name{dataset_num}]);
            datasets{dataset_num} = dataset;
            clear dataset
        end

    end

end
fprintf('done\n');



%% STEP 2: ???
fprintf('\n________________\nSTEP 2\n???\n');
dataset_num = 0;
% dataset_name = cell(1, length(animal_list)*length(area_list));
% datasets = cell(1, length(animal_list)*length(area_list));

for animal_num=1:length(animal_list)

    for area_num=1:length(area_list)

        dataset_num = dataset_num+1;
        fprintf('%i/%i\n', dataset_num, length(animal_list)*length(area_list));
        [datasets{dataset_num},collapsed{dataset_num}] = align_spike_data(datasets{dataset_num}, settings);

    end

end
fprintf('done\n');



%% STEP 3: ???
fprintf('\n________________\nSTEP 3\n???\n');

for animal_num=1:length(animal_list)

    for area_num=1:length(area_list)
        [datasets{dataset_num},collapsed{dataset_num}] = analyze_saccadic_responses(datasets{dataset_num}, collapsed{dataset_num}, settings);
    end

end