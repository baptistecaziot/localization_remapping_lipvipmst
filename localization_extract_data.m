
% Function to extract localization data
% Input:
% - monkey: 'hans' or 'chip'
% - area: 'LIP', 'VIP', 'MST' or 'MT'
% Output:
% - dataset: structure containing various fields, see readme for description


function dataset = localization_extract_data(monkey, area, settings)
    
    arguments
       monkey {mustBeText}
       area {mustBeText}
       settings struct
%        settings.spikes_filter_width (1,1) {mustBePositive} = 5
    end

    if strcmp(monkey, 'hans')
        monkey_folder = 'RawDataHans';
        extension = 'ha';
    elseif strcmp(monkey, 'chip')
        monkey_folder = 'RawDataChip';
        extension = 'ch';
    else
        error('Unknown monkey name');
    end
    if strcmp(area, 'LIP')
        area_folder = 'LIP_All';
    elseif strcmp(area, 'VIP')
        area_folder = 'VIP';
    elseif strcmp(area, 'MST')
        area_folder = 'MST';
    elseif strcmp(area, 'MT')
        area_folder = 'MT';
    elseif strcmp(area, 'MT_MST')
        area_folder = 'MT_MST';
    else
        error('Unknown area name');
    end

    data_path = ['../RFUp_Parietal_Compressed/', monkey_folder, '/', area_folder];
    sessions_list = dir([data_path,'/*.',extension,'d']);
    if isempty(sessions_list)
        error('Could not find data in folder %s', data_path);
    end

    fprintf('Found %i recording sessions\n', length(sessions_list))
    fprintf('Extracting data .. %03i%%', 0);

    for ses=1:length(sessions_list)
        fprintf('\b\b\b\b%03i%%', ceil(100*ses/length(sessions_list)));
        
        clear tmpsession
        
        tmpsession.path = data_path;
        tmpsession.extension = extension;
        [~,tmpsession.file] = fileparts(sessions_list(ses).name);

        tmpsession = read_data(tmpsession);
        tmpsession = read_indy(tmpsession);
        
        % Find common trials, this is because all trials are recorded in
        % the "stim_data" matrix, but only successful trials were recorded
        % in the neural data "spike_data".
        current_trial = 1;
        while 1
            if tmpsession.descr_mat(current_trial,2)~=tmpsession.stim_data(current_trial,3)
                nextind = 1;
                while 1
                    if ((current_trial+nextind)>size(tmpsession.descr_mat,1)) || ((current_trial+nextind)>size(tmpsession.stim_data,1))
                        tmpsession.analog_data(current_trial:end,:,:) = [];
                        tmpsession.spike_data(current_trial:end,:,:) = [];
                        tmpsession.events_time(current_trial:end,:) = [];
                        tmpsession.global_params(current_trial:end,:) = [];
                        tmpsession.descr_mat(current_trial:end,:) = [];
                        tmpsession.stim_data(current_trial:end,:) = [];
                        break;
                    end
                    if tmpsession.descr_mat(current_trial+nextind,2)==tmpsession.stim_data(current_trial,3)
                        tmpsession.analog_data(current_trial,:,:) = [];
                        tmpsession.spike_data(current_trial,:,:) = [];
                        tmpsession.events_time(current_trial,:) = [];
                        tmpsession.global_params(current_trial,:) = [];
                        tmpsession.descr_mat(current_trial,:) = [];
                        break;
                    elseif tmpsession.descr_mat(current_trial,2)==tmpsession.stim_data(current_trial+nextind,3)
                        tmpsession.stim_data(current_trial,:) = [];
                        break;
                    end
                    nextind = nextind+1;
                end
            else
                current_trial = current_trial+1;
            end
            
            if (current_trial>=size(tmpsession.descr_mat,1)) || (current_trial>=size(tmpsession.stim_data,1))
                break
            end
        end
        
        if size(tmpsession.stim_data,1)>size(tmpsession.descr_mat,1)
            tmpsession.stim_data(size(tmpsession.descr_mat,1)+1:end,:) = [];
        end
        tmpsession.trials_number = size(tmpsession.stim_data,1);
        
        if size(tmpsession.stim_data,1)~=size(tmpsession.descr_mat,1)
            error('Wrong trial number')
        end
        
        
        % Filter data
        filter_std = settings.spikes_filter_width;
        kernel = normpdf(-5*filter_std:+5*filter_std, 0, filter_std);
        for cc=1:7
            tmpsession.spike_data_filtered(:,:,cc) = settings.sr_neu * conv2(tmpsession.spike_data(:,:,cc), kernel, 'same');
        end

        dataset(ses) = tmpsession;
    end
    
    fprintf('\nAll done\n');

end



