
function [dataset,collapsed] = analyze_saccadic_responses(dataset, collapsed, settings)
%     
%     % -------------------------------------------------------------------------
%     % Script: analyze_stimulus_saccade_timing.m
%     %
%     % Description:
%     % This script analyzes the temporal alignment between visual stimulus onset,
%     % stimulus offset, and saccade onset across multiple experimental sessions.
%     % It compares corrected stimulus timings (from perisaccadic analysis) with
%     % original uncorrected values, visualizes distributions, and highlights
%     % potential misalignments or inconsistencies in timing.
%     %
%     % Key Features:
%     % - Computes and plots histograms of stimulus and saccade events
%     % - Visualizes trial-by-trial timing differences
%     % - Compares original vs corrected stimulus onset distributions
%     % - Performs session-level checks for missing or invalid data
%     %
%     % Requirements:
%     % - spike_data_info structure with perisaccadic_on_off field
%     % - dataset structure with stim_data field
%     % - condition_list and delay parameters
%     % -------------------------------------------------------------------------
%     
%     
%     % ---------------------------------------------------------------
%     % Initial Setup: Parameters and Data Structures
%     % ---------------------------------------------------------------
%     n_sessions   = spike_data_info.n_sessions;
%     n_conditions = spike_data_info.n_conditions;
%     delay        = spike_data_info.delay_ms; 
%     
%     % Precomputed perisaccadic timing data
%     perisaccadic_on_off = spike_data_info.perisaccadic_on_off;
%     
%     % ---------------------------------------------------------------
%     % Histogram: Absolute Difference Between Stimulus and Saccade Onset
%     % ---------------------------------------------------------------
%     figure;
%     all_times = [];
%     
%     for s = 1:n_sessions
%         if ~isempty(perisaccadic_on_off{s}) && istable(perisaccadic_on_off{s})
%             T = perisaccadic_on_off{s};
%             all_times = [all_times; T.stimulus_saccade_diff_ms];  % time difference in ms
%         end
%     end
%     
%     bin_width = 5;
%     histogram(all_times, 'BinWidth', bin_width); 
%     xlabel('abs(Stimulus Onset - Saccade Onset) (ms)');
%     ylabel('Count');
%     title('Distribution of abs(Stimulus Onset - Saccade Onset)');
%     grid on;
%     
%     % ---------------------------------------------------------------
%     % Histogram: Stimulus Onset, Offset, and Saccade Onset
%     % ---------------------------------------------------------------
%     figure;
%     hold on;
%     
%     all_onsets   = [];
%     all_offsets  = [];
%     all_saccades = [];
%     
%     for s = 1:n_sessions
%         if ~isempty(perisaccadic_on_off{s}) && istable(perisaccadic_on_off{s})
%             T = perisaccadic_on_off{s};
%             all_onsets   = [all_onsets; T.stimulus_onset_ms];
%             all_offsets  = [all_offsets; T.stimulus_offset_ms];
%             all_saccades = [all_saccades; T.saccade_onset_ms];
%         end
%     end
%     
%     % Define bin edges for consistent histogram binning
%     bin_width = 5;
%     min_edge = min([all_onsets; all_offsets; all_saccades]);
%     max_edge = max([all_onsets; all_offsets; all_saccades]);
%     bin_edges = min_edge:bin_width:max_edge;
%     
%     % Plot histograms with distinct colors
%     histogram(all_onsets,  bin_edges, 'FaceColor', [0 1 0],   'FaceAlpha', 0.6, 'DisplayName', 'Stimulus Onset');
%     histogram(all_offsets, bin_edges, 'FaceColor', [0 0.4 0], 'FaceAlpha', 0.6, 'DisplayName', 'Stimulus Offset');
%     histogram(all_saccades, bin_edges, 'FaceColor', [1 0 0],  'FaceAlpha', 0.6, 'DisplayName', 'Saccade Onset');
%     
%     xlabel('Time (ms)');
%     ylabel('Count');
%     title('Distribution of Stimulus and Saccade Events');
%     legend;
%     grid on;
%     
%     % ---------------------------------------------------------------
%     % Display Min/Max Ranges for Each Event Type
%     % ---------------------------------------------------------------
%     disp([min(all_onsets), max(all_onsets)]);
%     disp([min(all_offsets), max(all_offsets)]);
%     disp([min(all_saccades), max(all_saccades)]);
%     
%     % Per-session check for onset and saccade timing
%     for s = 1:n_sessions
%         T = perisaccadic_on_off{s};
%         
%         if isempty(T) || ~istable(T)
%             fprintf('Session %d: empty or not valid\n', s);
%             continue;
%         end
%     
%         disp([s, min(T.stimulus_onset_ms), min(T.saccade_onset_ms)]);
%     end
%     
%     % ---------------------------------------------------------------
%     % Scatter Plot: Stimulus-Saccade Timing Per Trial
%     % ---------------------------------------------------------------
%     figure;
%     hold on;
%     for s = 1:n_sessions
%         if ~isempty(perisaccadic_on_off{s}) && istable(perisaccadic_on_off{s})
%             T = perisaccadic_on_off{s};
%             scatter(1:length(T.stimulus_saccade_diff_ms), T.stimulus_saccade_diff_ms, 10, 'filled');
%         end
%     end
%     
%     xlabel('Trial Index');
%     ylabel('Stimulus Onset - Saccade Onset (ms)');
%     title('Stimulus-Saccade Timing Across Trials');
%     grid on;
%     
%     %%
%     % ---------------------------------------------------------------
%     % Extract Original (Uncorrected) Stimulus Onsets
%     % ---------------------------------------------------------------
%     onsets_original_all = cell(n_sessions,1);
%     
%     for s = 1:n_sessions
%         session = dataset(s);
%         onsets_stimulus_2_session = [];
%     
%         for c = 1:n_conditions
%             cond_id = condition_list(c);
%             trials = find(session.stim_data(:,3) == cond_id);
%     
%             for t = trials'  % loop over trial indices
%                 onset_stimulus_2 = session.stim_data(t,5) - delay;
%                 onsets_stimulus_2_session = [onsets_stimulus_2_session; onset_stimulus_2];
%             end
%         end
%         onsets_original_all{s} = onsets_stimulus_2_session;
%     end
%     
%     % Flatten all original onsets into one array
%     all_original_onsets = vertcat(onsets_original_all{:});
%     
%     % ---------------------------------------------------------------
%     % Recollect Perisaccadic Data for Comparison
%     % ---------------------------------------------------------------
%     all_onsets   = [];
%     all_offsets  = [];
%     all_saccades = [];
%     
%     for s = 1:n_sessions
%         if ~isempty(perisaccadic_on_off{s}) && istable(perisaccadic_on_off{s})
%             T = perisaccadic_on_off{s};
%             all_onsets   = [all_onsets; T.stimulus_onset_ms];
%             all_offsets  = [all_offsets; T.stimulus_offset_ms];
%             all_saccades = [all_saccades; T.saccade_onset_ms];
%         end
%     end
%     
%     % Define bin edges including original onsets
%     bin_width = 5;
%     min_edge = min([all_onsets; all_offsets; all_saccades; all_original_onsets]);
%     max_edge = max([all_onsets; all_offsets; all_saccades; all_original_onsets]);
%     bin_edges = min_edge:bin_width:max_edge;
%     
%     % ---------------------------------------------------------------
%     % Final Histogram: Original vs Corrected Stimulus Onsets
%     % ---------------------------------------------------------------
%     figure;
%     hold on;
%     
%     histogram(all_original_onsets, bin_edges, 'FaceColor', [0.2 0.8 1], 'FaceAlpha', 0.5, 'DisplayName', 'Original Stimulus Onset');
%     histogram(all_onsets,          bin_edges, 'FaceColor', [0 1 0],     'FaceAlpha', 0.6, 'DisplayName', 'Corrected Stimulus Onset');
%     histogram(all_offsets,         bin_edges, 'FaceColor', [0 0.4 0],   'FaceAlpha', 0.6, 'DisplayName', 'Stimulus Offset');
%     histogram(all_saccades,        bin_edges, 'FaceColor', [1 0 0],     'FaceAlpha', 0.6, 'DisplayName', 'Saccade Onset');
%     
%     xlabel('Time (ms)');
%     ylabel('Count');
%     title('Comparison of Original and Corrected Stimulus Onsets');
%     legend;
%     grid on;


    %% Check valid cells
    n_sessions = size(collapsed.peri.stim_sac_diff_mat,3);
    n_channels = size(collapsed.peri.stim_sac_diff_mat,2);
    n_trials = size(collapsed.peri.stim_sac_diff_mat,1);
    
    mean_fr = NaN(n_channels, n_sessions);
    valid_cells_mat = zeros(n_channels, n_sessions);
    valid_cells_count = 0;

    for ss=1:n_sessions
        for cc=1:n_channels
            
            mean_fr(cc,ss) = 1000*mean(collapsed.peri.spike_data_raw.saccade(:,:,:,cc,ss),'all','omitnan');  % mean across trials, time and condition
            valid_cells_mat(cc,ss) = mean_fr(cc,ss)>settings.min_FR;                                    % if it is higher than threshold
            valid_cells_count = valid_cells_count+1;

        end
    end
    
    fprintf('Valid cells: %i\n', valid_cells_count);
    
    collapsed.mean_fr = mean_fr;
    collapsed.valid_cells_mat = valid_cells_mat;
    
    
    %% Compute saccade-related activity, baseline etc.
    baseline_fr = NaN(n_channels, n_sessions);
    presaccade_fr = NaN(n_channels, n_sessions);
    saccade_fr = NaN(n_channels, n_sessions);
    postsaccade_fr = NaN(n_channels, n_sessions);

    n_figs = ceil(valid_cells_count / (settings.plot.n_rows*settings.plot.n_cols));
    fh_saccadeactivity = NaN(n_figs,1);
    for ff=1:n_figs
        fh_saccadeactivity(ff) = figure('Name', sprintf('Saccade-related activity (%i/%i)',ff,n_figs));
    end
    curr_cell = 0;
    for ss=1:n_sessions
        for cc=1:n_channels
            if valid_cells_mat(cc,ss)
                curr_cell = curr_cell+1;
                curr_fig = 1+fix((curr_cell-1)/(settings.plot.n_rows*settings.plot.n_cols));
                curr_ind = curr_cell-(curr_fig-1)*(settings.plot.n_rows*settings.plot.n_cols);
%                 curr_row = 1+fix(curr_ind/settings.plot.n_cols);
%                 curr_col = curr_ind-(curr_row-1)*settings.plot.n_cols;

                figure(fh_saccadeactivity(curr_fig))
                subplot(settings.plot.n_rows,settings.plot.n_cols,curr_ind)

                plot(settings.time_vec, mean(collapsed.sac.spike_data_filtered.saccade(:,:,1,cc,ss),1,'omitnan'));
                axis([-500,+500,0,200])

                % Baseline is relative to target jump
                baseline_fr(cc,ss) = 1000*mean(collapsed.sac.spike_data_raw.jump(:,(settings.time_vec>=settings.baseline_window(1))&(settings.time_vec<=settings.baseline_window(2)),1,cc,ss),'all','omitnan');
                % Rest relative to saccade
                presaccade_fr(cc,ss) = 1000*mean(collapsed.sac.spike_data_raw.saccade(:,(settings.time_vec>=settings.pre_saccadic_window(1))&(settings.time_vec<=settings.pre_saccadic_window(2)),1,cc,ss),'all','omitnan');
                saccade_fr(cc,ss) = 1000*mean(collapsed.sac.spike_data_raw.saccade(:,(settings.time_vec>=settings.saccadic_window(1))&(settings.time_vec<=settings.saccadic_window(2)),1,cc,ss),'all','omitnan');
                postsaccade_fr(cc,ss) = 1000*mean(collapsed.sac.spike_data_raw.saccade(:,(settings.time_vec>=settings.post_saccadic_window(1))&(settings.time_vec<=settings.post_saccadic_window(2)),1,cc,ss),'all','omitnan');

            end
        end
    end

    figure('Name','Windowed FR')
    subplot(2,3,1)
    plot([0,settings.fr_max],[0,settings.fr_max],'k--',...
        baseline_fr,saccade_fr,'ko')
    axis([0,settings.fr_max,0,settings.fr_max])
    xlabel('Baseline (Hz)')
    ylabel('Saccade (Hz)')

    subplot(2,3,2)
    plot([0,settings.fr_max],[0,settings.fr_max],'k--',...
        saccade_fr,presaccade_fr,'ko')
    axis([0,settings.fr_max,0,settings.fr_max])
    xlabel('Saccade (Hz)')
    ylabel('Pre-saccadic (Hz)')

    subplot(2,3,3)
    plot([0,settings.fr_max],[0,settings.fr_max],'k--',...
        saccade_fr,postsaccade_fr,'ko')
    axis([0,settings.fr_max,0,settings.fr_max])
    xlabel('Saccade (Hz)')
    ylabel('Post-saccadic (Hz)')

    subplot(2,3,4)
    hist(saccade_fr(:)-baseline_fr(:));

    subplot(2,3,5)
    hist(presaccade_fr(:)-saccade_fr(:));

    subplot(2,3,6)
    hist(postsaccade_fr(:)-saccade_fr(:));

    collapsed.baseline_fr = baseline_fr;
    collapsed.presaccade_fr = presaccade_fr;
    collapsed.saccade_fr = saccade_fr;
    collapsed.postsaccade_fr = postsaccade_fr;
    
end



