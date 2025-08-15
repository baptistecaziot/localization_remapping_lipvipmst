
function dataset = localization_eyeposition(dataset, options)
    
    % Options:
    % - eye (1, 2 or 3): compute eye position based of eye 1, 2 or both (3)
    % - threshold (>0): saccade speed threshold
    % - filter_width (int>0): width of convolution filter
    % - saccade_window_min: time-point to start clipping the saccades
    % - saccade_window_max: time-point to stop clipping the saccades
    
    arguments
       dataset struct
       options.eye (1,1) {mustBeMember(options.eye,[1,2,3])} = 3
       options.threshold (1,1) {mustBePositive} = 50
       options.filter_width (1,1) {mustBePositive,mustBeInteger} = 2
       options.saccade_window_min (1,1) = -100
       options.saccade_window_max (1,1) = +200
    end
    
    n_cells = length(dataset);
    n_col = 4;
    
    filter_vec = normpdf(-5*options.filter_width:+5*options.filter_width,0,options.filter_width);

    for cc=1:n_cells
        
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
        
        % Now find time-points where eye speed is higher than threshold and
        % after target onset
        indices_mat = repmat(1:size(dataset(cc).saccade.eye_speed,2),size(dataset(cc).saccade.eye_speed,1),1);
        fixation_onset_mat = repmat(dataset(cc).events_time(:,2),1,size(dataset(cc).saccade.eye_speed,2));
        
        saccade_indices = indices_mat.*((indices_mat>fixation_onset_mat) & (dataset(cc).saccade.eye_speed>options.threshold));
        saccade_indices(saccade_indices==0) = Inf;
        dataset(cc).saccade.saccade_onsets = min(saccade_indices, [] , 2);
        dataset(cc).saccade.saccade_latency = dataset(cc).saccade.saccade_onsets - dataset(cc).events_time(:,2);
        
        % Then we clip eye-position and eye-speed based on saccade onset,
        % for some reason matrix indexing doesn't work here so we'll do it
        % with a loop
        clipping_size = options.saccade_window_max-options.saccade_window_min+1;
        dataset(cc).saccade.saccade_pos = NaN(size(dataset(cc).saccade.saccade_onsets,1), clipping_size, 2);
        dataset(cc).saccade.saccade_speed = NaN(size(dataset(cc).saccade.saccade_onsets,1), clipping_size);
        for tt=1:size(dataset(cc).saccade.saccade_onsets,1)
            dataset(cc).saccade.saccade_pos(tt,:,:) = dataset(cc).saccade.eye_position(tt, dataset(cc).saccade.saccade_onsets(tt)+ (options.saccade_window_min:options.saccade_window_max), :);
            dataset(cc).saccade.saccade_speed(tt,:) = dataset(cc).saccade.eye_speed(tt, dataset(cc).saccade.saccade_onsets(tt)+ (options.saccade_window_min:options.saccade_window_max));
        end

        % Compute some other metrics
        dataset(cc).saccade.saccade_pos_start = squeeze(dataset(cc).saccade.saccade_pos(:,1,:));
        dataset(cc).saccade.saccade_pos_end = squeeze(dataset(cc).saccade.saccade_pos(:,end,:));
        dataset(cc).saccade.saccade_peak_velocity = max(dataset(cc).saccade.saccade_speed,[],2);

    end

end

