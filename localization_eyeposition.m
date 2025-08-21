
function dataset = localization_eyeposition(dataset, settings)
    
    % Options:
    % - eye (1, 2 or 3): compute eye position based of eye 1, 2 or both (3)
    % - threshold (>0): saccade speed threshold
    % - filter_width (int>0): width of convolution filter
    % - saccade_window_min: time-point to start clipping the saccades
    % - saccade_window_max: time-point to stop clipping the saccades
    
    arguments
       dataset (1,:) struct
       settings struct
    end

    options = settings.eyetracking;

    if ~isfield(options,'eye'); options.eye = 3; end
    if ~isfield(options,'speed_threshold_low_degsec'); options.speed_threshold_low_degsec = 50; end
    if ~isfield(options,'speed_threshold_high_degsec'); options.speed_threshold_high_degsec = 200; end
    if ~isfield(options,'filter_width'); options.filter_width = 2; end
    if ~isfield(options,'saccade_window_min'); options.saccade_window_min = -100; end
    if ~isfield(options,'saccade_window_max'); options.saccade_window_max = +200; end
    if ~isfield(options,'position_window_ms'); options.position_window_ms = 50; end

    n_cells = length(dataset);
    
    filter_vec = normpdf(-5*options.filter_width:+5*options.filter_width,0,options.filter_width);

    for cc=1:n_cells
        
        n_trials = size(dataset(cc).analog_data,1);

        % First convole eye-position data witha gaussian filter
        tmp = NaN(size(dataset(cc).analog_data));
        for ee=1:4
            tmp(:,:,ee) = conv2(dataset(cc).analog_data(:,:,ee), filter_vec, 'same');
        end
        
        % Compute mean of 2 eyes, or get data from only 1 eye
        if options.eye==3
            dataset(cc).saccade.eye_position(:,:,1) = mean(tmp(:,:,[1,3]),3);
            dataset(cc).saccade.eye_position(:,:,2) = mean(tmp(:,:,[2,4]),3);
        else
            dataset(cc).saccade.eye_position(:,:,1) = tmp(:,:,2*(options.eye-1)+1);
            dataset(cc).saccade.eye_position(:,:,2) = tmp(:,:,2*(options.eye-1)+2);
        end

        % Compute eye speed (diff of eye-position * 1000)
        dataset(cc).saccade.eye_speed = [zeros(size(tmp,1),1),1000*sqrt(diff(dataset(cc).saccade.eye_position(:,:,1),[],2).^2 + diff(dataset(cc).saccade.eye_position(:,:,1),[],2).^2)];
        
%         % Now find time-points after target onset where eye speed is higher than threshold
%         indices_mat = repmat(1:size(dataset(cc).saccade.eye_speed,2),size(dataset(cc).saccade.eye_speed,1),1);
%         fixation_onset_mat = repmat(dataset(cc).events_time(:,2),1,size(dataset(cc).saccade.eye_speed,2));
        
        
        % Compute saccade onset, offsets etc.
        clipping_size = options.saccade_window_max-options.saccade_window_min+1;
        time_vec = 1:size(dataset(cc).saccade.eye_speed,2);
        dataset(cc).saccade.bad_saccade = ones(n_trials,1);
        dataset(cc).saccade.saccade_onsets = NaN(n_trials,1);
        dataset(cc).saccade.saccade_offsets = NaN(n_trials,1);
        dataset(cc).saccade.saccade_latency = NaN(n_trials,1);
        dataset(cc).saccade.saccade_pos = NaN(n_trials, clipping_size, 2);
        dataset(cc).saccade.saccade_speed = NaN(n_trials, clipping_size);
        for tt=1:size(dataset(cc).saccade.saccade_onsets,1)
            try
                % Find saccade onset
                crossing_time = find((time_vec>dataset(cc).events_time(tt,2)) & (dataset(cc).saccade.eye_speed(tt,:)>options.speed_threshold_high_degsec), 1);
                dataset(cc).saccade.saccade_onsets(tt) = find(dataset(cc).saccade.eye_speed(tt,1:crossing_time)<options.speed_threshold_low_degsec,1,'last');
                dataset(cc).saccade.saccade_offsets(tt) = crossing_time + find(dataset(cc).saccade.eye_speed(tt,crossing_time:end)<options.speed_threshold_low_degsec,1) - 1;
                dataset(cc).saccade.saccade_latency = dataset(cc).saccade.saccade_onsets(tt)-dataset(cc).events_time(tt,2);
                
                dataset(cc).saccade.saccade_pos(tt,:,:) = dataset(cc).saccade.eye_position(tt, dataset(cc).saccade.saccade_onsets(tt) + (options.saccade_window_min:options.saccade_window_max), :);
                dataset(cc).saccade.saccade_speed(tt,:) = dataset(cc).saccade.eye_speed(tt, dataset(cc).saccade.saccade_onsets(tt) + (options.saccade_window_min:options.saccade_window_max));
                dataset(cc).saccade.bad_saccade(tt) = 0;
            catch
                if settings.verbose
                    fprintf('Failed to parse saccade\n');
                end
            end
        end
        
        % Compute some other metrics
        dataset(cc).saccade.saccade_pos_start = squeeze(median(dataset(cc).saccade.saccade_pos(:,1:options.position_window_ms,:),2));
        dataset(cc).saccade.saccade_pos_end = squeeze(median(dataset(cc).saccade.saccade_pos(:,end-options.position_window_ms:end,:),2));
        dataset(cc).saccade.saccade_peak_velocity = max(dataset(cc).saccade.saccade_speed,[],2);

    end

end

