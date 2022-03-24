function reliability_sampling_simulations(dat_in, cond_in, group_in)

% Currently setup to run through a SLURM-scheduled supercomputing cluster.
% To run normally, set SAMPLE to be the timepoint of interest and itterate
% through timepoints if necessary. Be aware this will take a long time if
% you are working off a single core and running serially.

% ex: {'hfa_pwr'}, {'targ_cor', 'dist_cor'; 'targ_cor_p1', 'dist_cor_p1'}, {0; 1:5; 6:10; 11:15}

summary.do_fast = 0;
itt_f = 1;

summary.dir_in = %%%%%%%%%%%%%%%%%%%%%%%%;
summary.dir_out = %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

summary.groupings = group_in;

summary.counts = [1:250];
summary.boots = 1000;
summary.crit_percent = 95;
summary.smooth_win = 2;

summary.d_type = dat_in;
summary.align_point = 'aot';

summary.conditions = cond_in;

boots_ref = summary.boots;
counts_ref = summary.counts;
crit_ref = summary.crit_percent;
groupings = summary.groupings;
set_size = 6;
smooth_win = summary.smooth_win;

SAMPLE = str2num(getenv('SLURM_ARRAY_TASK_ID'));
SAMPLE = SAMPLE + summary.smooth_win;

for I_D = 1 : size(summary.d_type, 2)
    for i_sw = 1:itt_f

        if itt_f > 1
            sessions_ind = find(fast_sessions==(i_sw-1));
        end

        for i_cp = 1 : size(summary.conditions, 1)

            load([summary.dir_in summary.d_type{I_D} '_' ...
                summary.conditions{i_cp,1} '_' summary.conditions{i_cp, 2} '_ext_v4.mat'], ...
                'col_pop_dist', 'col_pop_targ', 'depth_pop_dist', 'depth_pop_targ', ...
                'reli_dat_pop_dist', 'reli_dat_pop_targ');

            reli_dat_pop_targ = nanmean(reli_dat_pop_targ(:,SAMPLE-smooth_win:SAMPLE+smooth_win),2);
            reli_dat_pop_dist = nanmean(reli_dat_pop_dist(:,SAMPLE-smooth_win:SAMPLE+smooth_win),2);

            targ_nnan = ~isnan(reli_dat_pop_targ);
            dist_nnan = ~isnan(reli_dat_pop_dist);

            reli_dat_pop_targ = reli_dat_pop_targ(targ_nnan);
            reli_dat_pop_dist = reli_dat_pop_dist(dist_nnan);
            depth_pop_targ = depth_pop_targ(targ_nnan);
            depth_pop_dist = depth_pop_dist(dist_nnan);
            col_pop_targ = col_pop_targ(targ_nnan);
            col_pop_dist = col_pop_dist(dist_nnan);

            clear targ_nnan dist_nnan

            if itt_f > 1

                disp(numel(depth_pop_dist))
                depth_pop_dist = depth_pop_dist(ismember(col_pop_dist,sessions_ind));
                depth_pop_targ = depth_pop_targ(ismember(col_pop_targ,sessions_ind));
                reli_dat_pop_dist = reli_dat_pop_dist(ismember(col_pop_dist,sessions_ind));
                reli_dat_pop_targ = reli_dat_pop_targ(ismember(col_pop_targ,sessions_ind));
                col_pop_dist = col_pop_dist(ismember(col_pop_dist,sessions_ind));
                col_pop_targ = col_pop_targ(ismember(col_pop_targ,sessions_ind));
                disp(numel(depth_pop_dist))

            end

            output_targ_choice_nrs = nan(numel(counts_ref), size(groupings, 1));
            output_dist_choice_nrs = nan(numel(counts_ref), size(groupings, 1));
            output_crit_percent_nrs = nan(numel(counts_ref), size(groupings, 1));
            output_hit_crit_nrs = nan(numel(counts_ref), size(groupings, 1));

            for grp = 1 : size(groupings, 1)
                for counts = 1 : numel(counts_ref)

                    if mod(counts_ref(counts),size(groupings{grp},1)) ~= 0
                        continue
                    end

                    targ_choice_nrs = 0;
                    dist_choice_nrs = 0;

                    for boots = 1 : boots_ref

                        rsamps_nrs = [];

                        for n_grp = 1 : size(groupings{grp},1)

                            if ~isempty(find(depth_pop_targ >= groupings{grp}(n_grp, 1) & ...
                                    depth_pop_targ <= groupings{grp}(n_grp, end))) & ...
                                    ~isempty(find(depth_pop_dist >= groupings{grp}(n_grp, 1) & ...
                                    depth_pop_dist <= groupings{grp}(n_grp, end)))

                                temp_targ_grp = find(depth_pop_targ >= groupings{grp}(n_grp, 1) & ...
                                    depth_pop_targ <= groupings{grp}(n_grp, end));
                                temp_dist_grp = find(depth_pop_dist >= groupings{grp}(n_grp, 1) & ...
                                    depth_pop_dist <= groupings{grp}(n_grp, end));

                                if numel(temp_targ_grp) >= counts_ref(counts)/size(groupings{grp},1)

                                    temp_samps_nrs = [];

                                    for stim_row = 1 : set_size
                                        not_full = true;
                                        while not_full
                                            if stim_row == 1
                                                temp_samps_nrs(stim_row, :) = randsample(numel(temp_targ_grp), ...
                                                    counts_ref(counts)/size(groupings{grp},1), false);
                                                temp_samps_nrs(stim_row, :) = temp_targ_grp(temp_samps_nrs(stim_row, :));
                                            else
                                                temp_samps_nrs(stim_row, :) = randsample(numel(temp_dist_grp), ...
                                                    counts_ref(counts)/size(groupings{grp},1), false);
                                                temp_samps_nrs(stim_row, :) = temp_dist_grp(temp_samps_nrs(stim_row, :));
                                            end
                                            if sum(isnan(temp_samps_nrs(stim_row,:))) == 0
                                                not_full = false;
                                            end
                                        end
                                    end

                                    rsamps_nrs = cat(2, rsamps_nrs, temp_samps_nrs);

                                end
                            end
                        end

                        if ~isempty(rsamps_nrs)
                            stim_sum_nrs = nan(1,set_size);
                            for stim_itt = 1 : size(rsamps_nrs, 1)
                                if stim_itt == 1
                                    stim_sum_nrs(stim_itt) = ...
                                        sum(reli_dat_pop_targ(rsamps_nrs(stim_itt,:)));
                                else
                                    stim_sum_nrs(stim_itt) = ...
                                        sum(reli_dat_pop_dist(rsamps_nrs(stim_itt,:)));
                                end
                            end

                            if strcmp(summary.d_type, 'csd')
                                [~, maxind] =  min(stim_sum_nrs);
                            else
                                [~, maxind] =  max(stim_sum_nrs);
                            end

                            if maxind == 1; targ_choice_nrs = targ_choice_nrs + 1;
                            else; dist_choice_nrs = dist_choice_nrs + 1; end
                        end
                    end

                    if targ_choice_nrs | dist_choice_nrs

                        output_targ_choice_nrs(counts, grp) = ...
                            single(targ_choice_nrs);

                        output_dist_choice_nrs(counts, grp) = ...
                            single(dist_choice_nrs);

                        output_crit_percent_nrs(counts, grp) = ...
                            single(targ_choice_nrs / (targ_choice_nrs + dist_choice_nrs) * 100);

                        output_hit_crit_nrs(counts, grp) = ...
                            output_crit_percent_nrs(counts, grp) > crit_ref;

                    else

                        output_targ_choice_nrs(counts, grp) = ...
                            NaN;

                        output_dist_choice_nrs(counts, grp) = ...
                            NaN;

                        output_crit_percent_nrs(counts, grp) = ...
                            NaN;

                        output_hit_crit_nrs(counts, grp) = ...
                            NaN;

                    end
                end
            end

            if ~exist([summary.dir_out summary.d_type{I_D} '_' ...
                    summary.conditions{i_cp,1} '_' summary.conditions{i_cp, 2}], 'dir')

                mkdir([summary.dir_out summary.d_type{I_D} '_' ...
                    summary.conditions{i_cp,1} '_' summary.conditions{i_cp, 2}]);

            end

            save([summary.dir_out summary.d_type{I_D} '_' ...
                summary.conditions{i_cp,1} '_' summary.conditions{i_cp, 2} '_' ...
                'SAMPLE_' num2str(SAMPLE) '_resample.mat'], ...
                'summary','output_targ_choice_nrs', 'output_dist_choice_nrs', ...
                'output_crit_percent_nrs', 'output_hit_crit_nrs')

            clear col_pop_dist col_pop_targ depth_pop_dist depth_pop_targ reli_dat_pop_dist reli_dat_pop_targ

        end
    end
end