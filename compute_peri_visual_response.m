

function [dataset,collapsed] = compute_peri_visual_response(dataset, collapsed, settings)

    % Compute the visual response for periscaccadic stimuli by subtracting
    % the saccace-only response for each trial and compuing mean activity
    % within an activity window.

%     n_figs = ceil(valid_cells_count / (settings.plot.n_rows*settings.plot.n_cols));
%     fh_saccadeactivity = NaN(n_figs,1);
%     for ff=1:n_figs
%         fh_saccadeactivity(ff) = figure('Name', sprintf('Saccade-related activity (%i/%i)',ff,n_figs));
%     end
    
%     curr_cell = 0;
%     for ss=1:n_sessions
%         for cc=1:n_channels
%             if valid_cells_mat(cc,ss)
%                 curr_cell = curr_cell+1;
%                 curr_fig = 1+fix((curr_cell-1)/(settings.plot.n_rows*settings.plot.n_cols));
%                 curr_ind = curr_cell-(curr_fig-1)*(settings.plot.n_rows*settings.plot.n_cols);
% %                 curr_row = 1+fix(curr_ind/settings.plot.n_cols);
% %                 curr_col = curr_ind-(curr_row-1)*settings.plot.n_cols;
% 
%                 figure(fh_saccadeactivity(curr_fig))
%                 subplot(settings.plot.n_rows,settings.plot.n_cols,curr_ind)
% 
%                 plot(settings.time_vec, 1000*mean(collapsed.spike_data_filtered.saccade(:,:,7,cc,ss),1,'omitnan'));
% 
%                 % Baseline is relative to target jump
%                 baseline_fr(cc,ss) = 1000*mean(collapsed.spike_data_raw.jump(:,(settings.time_vec>=settings.baseline_window(1))&(settings.time_vec<=settings.baseline_window(2)),7,cc,ss),'all','omitnan');
%                 % Rest relative to saccade
%                 presaccade_fr(cc,ss) = 1000*mean(collapsed.spike_data_raw.saccade(:,(settings.time_vec>=settings.pre_saccadic_window(1))&(settings.time_vec<=settings.pre_saccadic_window(2)),7,cc,ss),'all','omitnan');
%                 saccade_fr(cc,ss) = 1000*mean(collapsed.spike_data_raw.saccade(:,(settings.time_vec>=settings.saccadic_window(1))&(settings.time_vec<=settings.saccadic_window(2)),7,cc,ss),'all','omitnan');
%                 postsaccade_fr(cc,ss) = 1000*mean(collapsed.spike_data_raw.saccade(:,(settings.time_vec>=settings.post_saccadic_window(1))&(settings.time_vec<=settings.post_saccadic_window(2)),7,cc,ss),'all','omitnan');
% 
%             end
%         end
%     end

end


