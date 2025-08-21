
function [dataset,collapsed] = align_spike_data(dataset, settings)

% =========================================================================
% Script: align_spike_data.m
%
% Description:
% This script performs preprocessing and alignment of spike data across
% multiple experimental sessions. It detects stimulus onsets, aligns spike
% trains to key behavioral events (stimulus, saccade, jump), bins the data
% into firing rates, and identifies perisaccadic trials based on timing
% thresholds. It also computes neuron activity masks based on firing rate
% thresholds and checks consistency across alignments.
% The output includes structured containers for further analysis and
% visualizations.
%
% Key Features:
% - Detects stimulus onset using prioritized channels (7 then 4)
% - Aligns spike data to stimulus, saccade, and jump events
% - Converts raw spike trains into binned firing rates (Hz)
% - Identifies perisaccadic trials based on stimulus-saccade timing
% - Excludes invalid trials with detailed logging
% - Computes active neuron masks per alignment
% - Verifies mask consistency across alignments
% - Summaries and metadata export
%
% Requirements:
% - 'dataset' structure with spike_data, saccade info, and event timings
% - Functions: compute_trial_indices, find_stim_onsets
%
% =========================================================================
    
    
    
    %% ------------------ User/session parameters --------------------------------------------
    % Sessions and neurons
    n_sessions = numel(dataset);
    n_neurons  = dataset(1).spike_channels_number;
    n_conditions = numel(settings.condition_list);
    
    % alignments to process
    alignments = {'saccade', 'stimulus', 'jump'};
    
    show_summary = true; % se false → nessuna stampa
    save_summary_to_file = true;
    summary_filename = 'pipeline_summary.txt';
    
    
    
    %% --- 0) Trial indices and sizes -----------------------------------------
    
    % Compute trial indices per (condition, session), and max per cell
    [trial_idx_cond, max_trials_per_cond_session] = compute_trial_indices(dataset, settings.condition_list);
    
    % Global maximum trials across all (cond x session) to preallocate [trial] dim in the 5D matrices 
    max_trial_cond = max(max_trials_per_cond_session(:));
    
    % Maximum trial length [samples] across sessions — to preallocate [time] dim
    max_trial_len_samples = max(arrayfun(@(x) size(x.spike_data, 2), dataset));
    
    
    %% --- 1) Stimulus detection (per session, trial) -------------------------
    % Detect stimulus onsets per trial using channel priority stim_channels.
    % The function also produces per-condition presence/absence for diagnostics.
    
    if settings.verbosity
        fprintf('Find stim onsets\n');
    end

    [trigger_data, trigger_report] = find_trigger_onsets(...
        dataset, settings.stim_channels, trial_idx_cond, settings);
    
    
    %% --- 2) Preallocate 5D containers ---------------------------------------
    % Shapes: [trial x time/bins x condition x neuron x session]
    % - spike_data_raw.(align): raw samples aligned (time dimension = samples)
    % - spike_data_win.(align): counts per bin (time dimension = bins)
    % - spike_data_binned.(align): firing rate (Hz), same shape as counts
    % - spike_data_valid.(align): boolean mask for valid bins
%     spike_data_raw    = struct();
%     spike_data_win    = struct();
%     spike_data_binned = struct();
%     spike_data_valid  = struct();
    spike_data_raw    = struct();
    spike_data_filtered = struct();
    
    peri_saccadic_mat =  NaN(max_trial_cond, n_conditions, n_sessions);

    for align_num = 1:numel(alignments)
        align_name = alignments{align_num};
%         spike_data_raw.(aling_name)    = NaN(max_trial_cond,    max_trial_len_samples,  n_conditions, n_neurons, n_sessions);
%         spike_data_win.(aling_name)    = NaN(max_trial_cond,    settings.n_bins,        n_conditions, n_neurons, n_sessions);
%         spike_data_binned.(aling_name) = NaN(max_trial_cond,    settings.n_bins,        n_conditions, n_neurons, n_sessions);
%         spike_data_valid.(aling_name)  = false(max_trial_cond,  settings.n_bins,        n_conditions, n_neurons, n_sessions);
        spike_data_raw.(align_name)      =  NaN(max_trial_cond, settings.n_bins, n_conditions, n_neurons, n_sessions);
        spike_data_filtered.(align_name) =  NaN(max_trial_cond, settings.n_bins, n_conditions, n_neurons, n_sessions);
    end
    
    % Containers for trial exclusions and reasons, per alignment/cond/session
    excluded_trial_indices_all = struct();
    excluded_trials_reasons_all = struct();
    for align_num = 1:numel(alignments)
        align_name = alignments{align_num};
        excluded_trial_indices_all.(align_name)  = cell(n_sessions, n_conditions);
        excluded_trials_reasons_all.(align_name) = cell(n_sessions, n_conditions);
    end
    
    
    %% --- 3) Align → bin → convert to Hz -------------------------------------
    
    if settings.verbosity
        fprintf('Clip neural data\n');
    end

    % --- Initialize container for perisaccadic timing info ---
    perisaccadic_on_off = cell(n_sessions, 1);
    
    for session_num=1:n_sessions
        session_dataset = dataset(session_num);
        session_trials_num = session_dataset.trials_number;
        trial_len = size(session_dataset.spike_data, 2);     % [samples per trial]
        n_neurons_session = size(session_dataset.spike_data, 3);
        
        if settings.verbosity>1
            fprintf('Session %d/%d : %d trials, %d samples per trial\n', session_num, n_sessions, session_trials_num, trial_len);
        end

        % Initialize session-specific container
        perisaccadic_on_off{session_num} = [];  
        
        for cond_num = 1:n_conditions
            cond_id = settings.condition_list(cond_num);
            cond_trials_idx = trial_idx_cond{cond_num, session_num};
            cond_trials_num = numel(cond_trials_idx);
            
            for align_num = 1:numel(alignments)
                align_name = alignments{align_num};
                excluded_trials_reasons = cell(cond_trials_num, 1);
                
                for tmp_trial_num = 1:cond_trials_num
                    tmp_trial_idx = cond_trials_idx(tmp_trial_num);
                    
                    try

                        % ---- 3.1 Determine alignment onset (in samples) ----
                        onset_samples = NaN;
                        
                        % Apply delay once to both onset and offset (in ms)
                        stim_onset_ms  = trigger_data(session_num).stim_idx{tmp_trial_idx}(1);
                        
                        trigg_offset = find(trigger_data(session_num).stim_idx{tmp_trial_idx}<trigger_data(session_num).stim_idx{tmp_trial_idx}(1)+settings.min_triggers_delay_ms,1);
                        if ~isempty(trigg_offset)
                            stim_offset_ms = trigg_offset;
                        else
                            stim_offset_ms = NaN;
                        end

                        % Calculate stimulus duration in milliseconds
                        stim_duration_ms = stim_offset_ms - stim_onset_ms;

                        % Get saccade onset
                        saccade_onset_ms = session_dataset.saccade.saccade_onsets(tmp_trial_idx);

                        % Compute absolute time difference between stimulus onset and saccade
                        delta_onset  = abs(stim_onset_ms - saccade_onset_ms);
                
                        % Mark trials as perisaccadic if either onset or offset is within ±200 ms
                        is_perisaccadic = (delta_onset <= settings.perisacc_thresh_ms);
                        
                        % Compute signed time difference (in ms) between stimulus offset and saccade onset
                        % Negative = stimulus before saccade, Positive = stimulus after saccade
                        stim_sacc_diff_ms = stim_onset_ms - saccade_onset_ms;
            
                        % Store: [stim_onset_ms, saccade_onset_ms, stim_offset_ms, stim_sacc_diff_ms, stim_duration_ms, is_perisaccadic, cond_id]
                        perisaccadic_on_off{session_num}(end+1, :) = [stim_onset_ms, saccade_onset_ms, stim_offset_ms, stim_sacc_diff_ms, stim_duration_ms, is_perisaccadic, cond_id, tmp_trial_idx];
    

                        switch align_name
                            case 'saccade'
    %                             if isfield(session_dataset, 'saccade') && isfield(session_dataset.saccade, 'saccade_onsets')
    %                                 if numel(session_dataset.saccade.saccade_onsets) >= tmp_trial_idx
                                        onset_samples = session_dataset.saccade.saccade_onsets(tmp_trial_idx);
    %                                 end
    %                             end
                            
                            case 'stimulus'
%                                 % Check if stimulus indices are valid
%                                 if tmp_trial_idx <= numel(trigger_data(session_num).stim_idx) && ...
%                                     ~isempty(trigger_data(session_num).stim_idx{tmp_trial_idx}) && ...
%                                     ~isnan(trigger_data(session_num).stim_idx{tmp_trial_idx}(1)) && ...
%                                     numel(trigger_data(session_num).stim_idx{tmp_trial_idx}) >= 2 && ...
%                                     ~isnan(trigger_data(session_num).stim_idx{tmp_trial_idx}(2))
                                
                                if isempty(trigger_data(session_num).stim_idx{tmp_trial_idx}) ||...
                                        isnan(trigger_data(session_num).stim_idx{tmp_trial_idx}(1)) ||...
                                        isnan(session_dataset.saccade.saccade_onsets(tmp_trial_idx))
%                                     fprintf('Missing trigger\n');
                                    stim_sacc_diff_ms = NaN;
                                    continue;

                                else

                                    onset_samples = stim_onset_ms;
                                end
        
                            case 'jump'
                                if isfield(session_dataset, 'events_time') && size(session_dataset.events_time, 1) >= tmp_trial_idx
                                    onset_samples = session_dataset.events_time(tmp_trial_idx, 2); % event column 2 = jump --> its written int 
                                end
                        end


%                         % ---- 3.2 Exclude if onset invalid ----
%                         excluded_trials = zeros(n_conditions, n_sessions);
%                         trial_reasons = {};
%                         if ~(isnumeric(onset_samples) && isfinite(onset_samples) && (onset_samples>=1))
%                             trial_reasons{end+1} = sprintf('no_onset_%s', alg);
%                             excluded_trials_reasons{tmp_trial_num} = trial_reasons;
%                             excluded_trials(cond_num, session_num) = excluded_trials(cond_num, session_num)+1;
%                             continue;
%                         end
                        
    
                        % ---- 3.3 Compute window indices (in samples) ----
                        window_start_samples = round(settings.time_window(1) * settings.sr_neu / 1000); % [samples]
                        window_end_samples   = round(settings.time_window(2) * settings.sr_neu / 1000); % [samples]
        
                        % If you adopt the safer delay handling, adjust here:
                        % start_idx = round(onset_samples + window_start_samples - (an=="stimulus")*delay_samples);
                        % end_idx   = round(onset_samples + window_end_samples   - 1 - (an=="stimulus")*delay_samples);
                        start_idx = round(onset_samples + window_start_samples);
                        end_idx   = round(onset_samples + window_end_samples - 1);
        
                        % Exclude only if window has NO overlap with the trial at all
                        if end_idx < 1 || start_idx > trial_len
                            trial_reasons{end+1} = sprintf('window_out_of_range_%s', align_name);
                            excluded_trials_reasons{tmp_trial_num} = trial_reasons;
                            excluded_trials(cond_num, session_num) = excluded_trials(cond_num, session_num)+1;
                            continue;
                        end
                        
    
                        % ---- 3.4 Clamp window to valid bounds and build aligned array ----
                        valid_start = max(1, start_idx);
                        valid_end   = min(trial_len, end_idx);
        
                        % Expected aligned window length [samples]
                        expected_len = window_end_samples - window_start_samples;
        
                        % Allocate aligned raw buffer (zeros → implicit padding where outside)
                        aligned_raw = NaN(expected_len, n_neurons_session);
                        aligned_filtered = NaN(expected_len, n_neurons_session);
        
                        % Compute insertion positions inside the aligned buffer
                        insert_start = valid_start - start_idx + 1;
                        insert_end   = insert_start + (valid_end - valid_start);
        
%                         % Safety: clamp insertion indices to buffer limits
%                         if insert_start < 1
%                             valid_start = valid_start + (1 - insert_start);
%                             insert_start = 1;
%                         end
%                         if insert_end > expected_len
%                             valid_end = valid_end - (insert_end - expected_len);
%                             insert_end = expected_len;
%                         end
        
%                         if valid_start <= valid_end && insert_start <= insert_end
                            % Extract raw segment [samples x neurons]
                            clipped_spikes = squeeze(session_dataset.spike_data(tmp_trial_idx, valid_start:valid_end, :));
                            clipped_spikes_filtered = squeeze(session_dataset.spike_data_filtered(tmp_trial_idx, valid_start:valid_end, :));
        
%                             % Ensure rows = time, columns = neurons
%                             if size(seg, 1) ~= (valid_end - valid_start + 1) && size(seg, 2) == (valid_end - valid_start + 1)
%                                 seg = seg';
%                             end
        
                            % Place segment into aligned buffer
                            aligned_raw(insert_start:insert_end, :) = clipped_spikes;
                            aligned_filtered(insert_start:insert_end, :) = clipped_spikes_filtered;
%                         end
        
%                         % Ensure aligned_raw matches preallocated time dimension
%                         if size(aligned_raw, 1) > max_trial_len_samples
%                             aligned_raw = aligned_raw(1:max_trial_len_samples, :);
%                         elseif size(aligned_raw, 1) < max_trial_len_samples
%                             tmp = zeros(max_trial_len_samples, n_neurons_session);
%                             tmp(1:size(aligned_raw, 1), :) = aligned_raw;
%                             aligned_raw = tmp;
%                         end
        
                        % Store aligned raw
                        dataset(session_num).neural.raw.(align_name)(tmp_trial_num, :, cond_num, 1:n_neurons_session) = aligned_raw;
                        dataset(session_num).neural.filtered.(align_name)(tmp_trial_num, :, cond_num, 1:n_neurons_session) = aligned_filtered;
        
%                         % ---- 3.5 Binning: counts and Hz ----
%                         samples_for_bins = settings.bin_size_samples * settings.n_bins; % required total samples (bin_size_samples)
%                         data_to_bin = aligned_raw(1:min(samples_for_bins, size(aligned_raw, 1)), :);
%         
%                         % Pad with zeros if needed (top-up to exact multiple)
%                         if size(data_to_bin, 1) < samples_for_bins
%                             tmp = zeros(samples_for_bins, n_neurons_session);
%                             tmp(1:size(data_to_bin, 1), :) = data_to_bin;
%                             data_to_bin = tmp;
%                         end
%         
%                         bin_counts = zeros(settings.n_bins, n_neurons_session);
%                         for n = 1:n_neurons_session
%                             vec = data_to_bin(:, n);
%                             % Safety assertion for reshape
%                             assert(numel(vec) == settings.bin_size_samples * settings.n_bins, 'Binning size mismatch.'); % (bin_size_samples)
%                             vec_reshaped = reshape(vec, settings.bin_size_samples, settings.n_bins); % (bin_size_samples)
%                             bin_counts(:, n) = sum(vec_reshaped, 1).';
%                         end
                        
%                         % Store counts and validity
%                         spike_data_win.(alg)(tmp_trial_num, :, cond_num, 1:n_neurons_session, session_num)   = bin_counts;
%                         spike_data_valid.(alg)(tmp_trial_num, :, cond_num, 1:n_neurons_session, session_num) = true;
%                         
%                         % Convert to Hz (spikes per second)
%                         bin_duration_s = settings.bin_size_ms / 1000; % [s]
%                         spike_data_binned.(alg)(tmp_trial_num, :, cond_num, 1:n_neurons_session, session_num) = bin_counts / bin_duration_s;
                        
                        peri_saccadic_mat(tmp_trial_num, cond_num, session_num) = stim_sacc_diff_ms;
                        spike_data_raw.(align_name)(tmp_trial_num, :, cond_num, 1:n_neurons_session, session_num) = aligned_raw;
                        spike_data_filtered.(align_name)(tmp_trial_num, :, cond_num, 1:n_neurons_session, session_num) = aligned_filtered;
                    
                    catch ME
                        fprintf(ME.message);
                        keyboard;
                    end % try

                end % trials
    
%                 % Save exclusion bookkeeping
%                 excluded_trials_reasons_all.(alg){session_num, cond_num} = excluded_trials_reasons;
%                 excluded_mask_local = cellfun(@(x) ~isempty(x), excluded_trials_reasons);
%                 excluded_trial_indices_all.(alg){session_num, cond_num} = cond_trials_idx(excluded_mask_local);
    
                % Optional: exclusion report
                % n_excluded = sum(excluded_mask_local);
                % if n_excluded > 0
                %     fprintf('Excluded %d/%d trials for %s, cond %d\n', n_excluded, n_trials_cond, an, cond_id);
                % end
            end % alignments
        end % conditions
    
%         % Convert each session's matrix to a labeled table
%        if ~isempty(perisaccadic_on_off{session_num})
%             perisaccadic_on_off{session_num} = array2table(perisaccadic_on_off{session_num}, ...
%                 'VariableNames', {'stimulus_onset_ms', 'saccade_onset_ms', 'stimulus_offset_ms', 'stimulus_saccade_diff_ms', 'stimulus_duration_ms', 'is_perisaccadic', 'condition_id', 'trial_idx'});
%        end
        
    end % sessions
    
    
%     %% --- 4) Active neurons per alignment ------------------------------------
%     % Compute an activity mask [session x neuron] per alignment, excluding stim channels.
%     active_neurons = struct();
%     for align_num = 1:numel(alignments)
%         alg = alignments{align_num};
%         mask_active = false(n_sessions, n_neurons);
%         if settings.verbosity>1
%             fprintf('Computing active neurons for %s alignment...\n', alg);
%         end
%     
%         for session_num = 1:n_sessions
%             for n = 1:n_neurons
%                 % Skip stim channels for FR metrics
%                 if ismember(n, settings.stim_channels)
%                     mask_active(session_num, n) = false;
%                     continue;
%                 end
%     
%                 % Aggregate FR across all valid trials/conds for this neuron/session
%                 FR_all_trials = [];
%                 for cond_num = 1:n_conditions
%                     trials_loc = trial_idx_cond{cond_num, session_num};
%                     n_trials_loc = numel(trials_loc);
%                     if n_trials_loc == 0, continue; end
%     
%                     data_loc  = squeeze(spike_data_binned.(alg)(1:n_trials_loc, :, cond_num, n, session_num)); % [trials x bins]
%                     valid_loc = squeeze(spike_data_valid.(alg)(1:n_trials_loc, :, cond_num, n, session_num));   % [trials x bins]
%                     valid_trials = any(valid_loc, 2);
%                     if any(valid_trials)
%                         FR_all_trials = [FR_all_trials; data_loc(valid_trials, :)]; %#ok<AGROW>
%                     end
%                 end
%                 
%                 if isempty(FR_all_trials), continue; end
%                 mean_FR = mean(FR_all_trials(:), 'omitnan');
%                 mask_active(session_num, n) = (~isnan(mean_FR)) && (mean_FR >= settings.min_FR);
%             end
%         end
%     
%         active_neurons.(alg) = mask_active;
%         if settings.verbosity>1
%             fprintf('  Active neurons: %d/%d\n', sum(mask_active(:)), numel(mask_active));
%         end
%     end
    
%     %% --- 5) Mask consistency across alignments ------------------------------
%     % If identical across alignments → use a single mask; otherwise keep struct.
%     align_names = fieldnames(active_neurons);
%     base_mask = active_neurons.(align_names{1});
%     is_identical_across_alignments = true;
%     for k = 2:numel(align_names)
%         if ~isequal(base_mask, active_neurons.(align_names{k}))
%             is_identical_across_alignments = false;
%             break;
%         end
%     end
%     
%     if is_identical_across_alignments
%         fprintf('Active neuron mask identical across alignments. Using single mask.\n');
%         final_active_neurons = base_mask;        % [n_sessions x n_neurons]
%     else
%         fprintf('Active neuron mask differs across alignments. Keeping per-alignment masks.\n');
%         final_active_neurons = active_neurons;   % struct of masks
%     end
%     
%     fprintf('\nPipeline finished.\n');
%     
% 
%     %% --- 6) Summaries (all or nothing) --------------------------------------
%     
%     
%     if show_summary
%         if save_summary_to_file
%             fid = fopen(summary_filename, 'w');
%             outfun = @(varargin) fprintf(fid, varargin{:});
%         else
%             outfun = @fprintf;
%         end
%     
%         outfun('\n== DATA SUMMARY ==\n');
%         % Matrix dimensions
%         for align_num = 1:numel(alignments)
%             alg = alignments{align_num};
%             sz_raw    = size(spike_data_raw.(alg));
%             sz_win    = size(spike_data_win.(alg));
%             sz_binned = size(spike_data_binned.(alg));
%             outfun('%-9s | raw %s | win %s | bin %s\n', alg, mat2str(sz_raw), mat2str(sz_win), mat2str(sz_binned));
%         end
%     
%         % Percentage of valid bins
%          outfun('\n== Percentage of valid bins ==\n');
%         for align_num = 1:numel(alignments)
%             alg = alignments{align_num};
%             valid_mask = spike_data_valid.(alg);
%             pct_valid = 100 * nnz(valid_mask) / numel(valid_mask);
%             outfun('%-9s | %5.1f%% valid\n', alg, pct_valid);
%         end
%     
%         % FR stats per condition
%          outfun('\n== Firing rate (Hz) by condition ==\n');
%         for align_num = 1:numel(alignments)
%             alg = alignments{align_num};
%             outfun('%s:\n', alg);
%             for cond_num = 1:n_conditions
%                 X = squeeze(spike_data_binned.(alg)(:, :, cond_num, :, :));
%                 X = X(:);
%                 X = X(~isnan(X));
%                 if isempty(X)
%                      outfun('  cond %d | no data\n', settings.condition_list(cond_num));
%                 else
%                      outfun('  cond %d | min %.1f, max %.1f, mean %.1f Hz\n', ...
%                         settings.condition_list(cond_num), min(X), max(X), mean(X));
%                 end
%             end
%         end
%     
%         % Stimulus presence summary (per session/condition)
%          outfun('\n== Stimulus summary per condition ==\n');
%         for session_num = 1:numel(trigger_data)
%              outfun('Session %d:\n', session_num);
%             for cond_num = 1:numel(trigger_data(session_num).per_condition)
%                 pc = trigger_data(session_num).per_condition(cond_num);
%                  outfun('  cond %d | trials: %d | stim yes: %d | stim no: %d\n', ...
%                     pc.cond_id, numel(pc.trials), pc.n_with_stim, pc.n_without);
%             end
%         end
%     
%         % Compact stimulus channel usage
%         outfun('\n== Stimulus channel usage (compact) ==\n');
%         outfun('Session | Ch%-3d  Ch%-3d  NoStim\n', settings.stim_channels(1), settings.stim_channels(2));
%         outfun('--------+------------------------\n');
%         for session_num = 1:numel(trigger_data)
%             src   = trigger_data(session_num).source_channel;
%             n_ch1 = sum(src == settings.stim_channels(1));
%             n_ch2 = sum(src == settings.stim_channels(2));
%             n_ns  = sum(isnan(src));
%             outfun('  %2d    | %4d   %4d   %4d\n', session_num, n_ch1, n_ch2, n_ns);
%         end
%     
%         if save_summary_to_file
%             fclose(fid);
%             fprintf('Summary saved to %s\n', summary_filename);
%         end
%     end
    
    
    
%     
%     %% ------------------ Metadata export (spike_data_info) -------------------
%     spike_data_info = struct();
%     
%     % Basic experiment info
%     spike_data_info.animal     = animal;
%     spike_data_info.area       = area;
%     
%     % Sessions & neurons
%     spike_data_info.n_sessions = n_sessions;
%     spike_data_info.n_neurons  = n_neurons;
%     
%     % Conditions
%     spike_data_info.condition_list = settings.condition_list;
%     spike_data_info.n_conditions   = n_conditions;
%     spike_data_info.stim_positions = stim_positions;
%     spike_data_info.max_trial_cond = max_trial_cond;
%     
%     % Sampling & time (units annotated)
%     spike_data_info.sr_eye           = sr_eye;                   % [Hz]
%     spike_data_info.sr_neu           = sr_neu;                   % [Hz]
%     spike_data_info.time_window_ms   = settings.time_window;              % [ms]
%     spike_data_info.time_vec_ms      = settings.time_vec;                 % [ms]
%     spike_data_info.delay_ms         = delay_ms;                 % [ms]
%     spike_data_info.delay_samples    = delay_samples;            % [samples]
%     spike_data_info.bin_size_ms      = settings.bin_size_ms;    % [ms]
%     spike_data_info.bin_size_samples = bin_size_samples;         % [samples]
%     spike_data_info.n_bins           = settings.n_bins;
%     
%     % Thresholds and alignments
%     spike_data_info.min_FR_Hz        = min_FR;                   % [Hz]
%     spike_data_info.alignments       = alignments;
%     spike_data_info.perisaccadic_on_off = perisaccadic_on_off;   % [samples]
%     
%     % Active neuron mask info
%     spike_data_info.is_identical_across_alignments = is_identical_across_alignments;
%     
%     % Plotting
%     spike_data_info.subplot_map = subplot_map;
%     spike_data_info.cond_keys   = keys(spike_data_info.subplot_map);
%     spike_data_info.colors      = lines(length(spike_data_info.cond_keys));
%     spike_data_info.color_spike = color_spike;
%     spike_data_info.color_stim  = color_stim;
%     spike_data_info.color_jump  = color_jump;
%     spike_data_info.color_sacc  = color_sacc;
%     spike_data_info.size_tick   = size_tick;
    
    
    
%     %%% validation_script.m
%     fprintf('\n== VALIDATION CHECKS ==\n');
%     for align_num = 1:numel(alignments)
%         alg = alignments{align_num};
%         raw_sz = size(spike_data_raw.(alg));
%         win_sz = size(spike_data_win.(alg));
%         bin_sz = size(spike_data_binned.(alg));
%     
%         fprintf('\nAlignment: %s\n', alg);
%         fprintf('  raw size: %s\n', mat2str(raw_sz));
%         fprintf('  win size: %s\n', mat2str(win_sz));
%         fprintf('  bin size: %s\n', mat2str(bin_sz));
%     
%         % Check NaN patterns: raw should have NaNs beyond trial count
%         nan_ratio_raw = mean(isnan(spike_data_raw.(alg)(:)));
%         nan_ratio_bin = mean(isnan(spike_data_binned.(alg)(:)));
%         fprintf('  NaN ratio raw: %.2f%%\n', nan_ratio_raw*100);
%         fprintf('  NaN ratio bin: %.2f%%\n', nan_ratio_bin*100);
%         
%         % Check alignment between valid mask and NaNs in FR
%         valid_mask = spike_data_valid.(alg);
%         mismatch = xor(valid_mask, ~isnan(spike_data_binned.(alg)));
%         fprintf('  Mismatch valid vs NaNs: %d elements\n', nnz(mismatch));
%     end
%     fprintf('\nValidation complete.\n');
    
    
    
%     %% ------------------ Clear workspace (safe) ------------------------------
%     % Keep only what downstream analysis will need.
%     clearvars -except ...
%         dataset trial_idx_cond ...
%         trigger_data stim_report spike_data_info...
%         spike_data_raw spike_data_win spike_data_binned spike_data_valid ...
%         excluded_trials_reasons_all excluded_trial_indices_all ...
%         active_neurons final_active_neurons ...

%     for align_num = 1:numel(alignments)
%         align_name = alignments{align_num};

    collapsed.peri_saccadic_mat = peri_saccadic_mat;
    collapsed.spike_data_raw = spike_data_raw;
    collapsed.spike_data_filtered = spike_data_filtered;

end



%% ------------------ Helper functions ------------------------------------
function [trial_idx_cond, max_trials_per_cond_session] = compute_trial_indices(dataset, condition_list)
    % compute_trial_indices
    % Returns:
    % - trial_idx_cond: cell [n_conditions x n_sessions], trial indices per condition/session
    % - max_trials_per_cond_session: [n_conditions x n_sessions] counts per cell
    n_sessions   = numel(dataset);
    n_conditions = numel(condition_list);
    trial_idx_cond = cell(n_conditions, n_sessions);
    max_trials_per_cond_session = zeros(n_conditions, n_sessions);
    
    for s = 1:n_sessions
        session = dataset(s);
        for c = 1:n_conditions
            cond_id = condition_list(c);
            trials = find(session.stim_data(:, 3) == cond_id); % assumes col 3 holds condition ID
            trial_idx_cond{c, s} = trials;
            max_trials_per_cond_session(c, s) = numel(trials);
        end
    end
end

function [trigger_data, trigger_report] = find_trigger_onsets(dataset, trigger_channels, trial_idx_cond, settings)
    % find_stim_onsets
    % Detect TTL-based trigger onsets per trial and compute per-condition presence stats.
    % Priority: trigger_channels(1) first, then trigger_channels(2).
    %
    % Outputs (per session s):
    % - trigger_data(s).stim_idx{t}: vector of sample indices where TTL == 1, or NaN if none
    % - trigger_data(s).no_stim_trials(t): logical, true if no TTL found
    % - trigger_data(s).source_channel(t): channel index that provided TTL (or NaN)
    % - trigger_data(s).n_valid: number of trials with detected TTL
    % - trigger_data(s).per_condition(c): struct with cond_id, trials, stim_present[], n_with_stim, n_without
    % - report: aggregate counts across sessions (total/with/without)
    n_sessions = numel(dataset);
    n_conditions = numel(settings.condition_list);
    
    trigger_data = struct();
    trigger_report = struct('n_total', 0, 'n_with_stim', 0, 'n_without', 0);
    
    for session_num = 1:n_sessions
        session = dataset(session_num);
        n_trials = session.trials_number;
        stim_idx = cell(n_trials, 1);
        no_trigger_trials = false(n_trials, 1);
        source_channel = NaN(n_trials, 1);
        
        if settings.verbosity>1
            fprintf('Searching for stimulus onsets in session %d/%d...\n', session_num, n_sessions);
        end

        for trial = 1:n_trials
            found = false;
            for ch_try = 1:numel(trigger_channels)
                ch = trigger_channels(ch_try);
                if ch > size(session.spike_data, 3)
                    continue; % channel absent
                end
                stim_vec = squeeze(session.spike_data(trial, :, ch));
                idx = find(stim_vec == 1); % binary TTL expected
                if ~isempty(idx)
                    stim_idx{trial} = idx;
                    source_channel(trial) = ch;
                    found = true;
                    break; % respect priority order
                end
            end
            if ~found
                stim_idx{trial} = NaN;
                no_trigger_trials(trial) = true;
            end
        end
    
        trigger_data(session_num).stim_idx       = stim_idx;
        trigger_data(session_num).no_stim_trials = no_trigger_trials;
        trigger_data(session_num).source_channel = source_channel;
        trigger_data(session_num).n_valid        = sum(~no_trigger_trials);
    
        % Per-condition presence/absence
        per_condition = struct();
        for c = 1:n_conditions
            trials_c = trial_idx_cond{c, session_num};
            stim_present = false(size(trials_c));
            for k = 1:numel(trials_c)
                trial = trials_c(k);
                stim_present(k) = ~isempty(stim_idx{trial}) && ~isnan(stim_idx{trial}(1));
            end
            per_condition(c).cond_id      = settings.condition_list(c);
            per_condition(c).trials       = trials_c;
            per_condition(c).stim_present = stim_present;
            per_condition(c).n_with_stim  = sum(stim_present);
            per_condition(c).n_without    = sum(~stim_present);
        end
        trigger_data(session_num).per_condition = per_condition;
        if settings.verbosity>1
            fprintf('  Found stimulus in %d/%d trials\n', trigger_data(session_num).n_valid, n_trials);
        end
    end
    
    % Aggregate report across sessions
    trigger_report.n_total     = sum(arrayfun(@(x) x.trials_number, dataset));
    trigger_report.n_with_stim = sum(arrayfun(@(sd) sd.n_valid, trigger_data));
    trigger_report.n_without   = trigger_report.n_total - trigger_report.n_with_stim;
    
    if settings.verbosity>1
        fprintf('Stimulus onset search completed: %d trials with stim, %d without.\n', ...
            trigger_report.n_with_stim, trigger_report.n_without);
    end
end

