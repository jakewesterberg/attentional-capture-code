function summary = compile_data()

% Jake Westerberg
% Vanderbilt University
% jakewesterberg@gmail.com

%% INPUT PARAMETERS
summary.dir_in = 'Z:\_BACKUPS\Xwings_backup\V4\VS\';
summary.dir_out = 'Z:\_DATA\MANUSCRIPTS\WesterbergEtAl_2022_NCOMMS_Feedforward\Revision_Data\';
summary.d_type = {'mua'};
summary.align_point = 'aot';

summary.conditions = { ...
    'targ_cor', 'dist_adj_cor'; ...
    'targ_cor', 'dist_opp_cor'; ...
    'targ_cor_p0', 'dist_cor_p0';  ...
    'targ_cor_p1', 'dist_cor_p1';  ...
    ...'targ_err', 'dist_err_sac' ...
    ...'targ_cor', 'dist_cor' ...
    };
    
%% ANALYSIS
good_files = [1 2	4	5	6	7	8	9	11	13	14	15	16	17	18	22	24	27	28	30	31	34	35	38	39	40	41	42	44	44];
good_probes = [1 1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	1	2];

[rec, rec_head] = H_REC();
dir_dir = H_DIRF(summary.dir_in);
summary.directories = dir_dir;

if strcmp(summary.align_point, 'aos')
    comp_pre = 152;     comp_post = 152;
else
    comp_pre = 105;     comp_post = 301;
end

for id = 1 : size(summary.d_type, 2)
    for cp = 1 : size(summary.conditions, 1)

        reli_dat_pop_targ = [];     reli_dat_pop_dist = [];
        depth_pop_targ = [];        depth_pop_dist = [];
        col_pop_targ = [];          col_pop_dist = [];
        monk_pop_targ = [];         monk_pop_dist = [];    
        rt_targ = [];               rt_dist = [];
        pr_targ = [];               pr_dist = [];
        rtrnk_targ = [];            rtrnk_dist = [];
        dist_from_targ = [];    
        
        p_ctr = 0;
        
        for i_dir = 1 : numel(dir_dir)
            
            analy.cur_dir = dir_dir{i_dir};
            
            idx_date  = strfind(analy.cur_dir, '-');
            analy.date  = str2double(analy.cur_dir(idx_date(end-1)+1:idx_date(end)-1));
            analy.time = str2double(analy.cur_dir(idx_date(end)+1:idx_date(end)+6));
            
            idx_slash  = strfind(analy.cur_dir, '\');
            analy.file_ident = analy.cur_dir(idx_slash(end)+1:end);
            
            analy.ind_prb = find( analy.date == [rec{:,rec_head.DATE}]);
            
            for i_prb = 1 : length(analy.ind_prb)
                
                ident_file = find(good_files(good_files == i_dir));
                
                if sum(ident_file) == 0
                    continue
                else
                    is_it_a_good_probe = sum(i_prb == good_probes(ident_file));
                end
                
                if strcmp(rec{ analy.ind_prb(i_prb), rec_head.PROCESS }, 'exclude') | ...
                        ~is_it_a_good_probe
                    continue
                else
                    p_ctr = p_ctr + 1;
                end
                
                disp([ 'COMPILING DATA: Probe ' num2str(p_ctr) ])
                
                analy.rf_loc = rec{ analy.ind_prb(i_prb), rec_head.RF };
                
                analy.L4_top = rec{ analy.ind_prb(i_prb), rec_head.TO4 };
                analy.L4_bot = rec{ analy.ind_prb(i_prb), rec_head.BO4 };
                analy.Cx_top = rec{ analy.ind_prb(i_prb), rec_head.TOCX };
                analy.Cx_bot = rec{ analy.ind_prb(i_prb), rec_head.BOCX };
                
                analy.probes = rec{ analy.ind_prb(i_prb), rec_head.TPROBE };
                analy.location = rec{ analy.ind_prb(i_prb), rec_head.LOCATION};           
                analy.monkey = rec{ analy.ind_prb(i_prb), rec_head.IDENT };
                
                if strcmp(rec{ analy.ind_prb(i_prb), rec_head.SORTDIR }, 'descending')
                    analy.ch_min = rec{ analy.ind_prb(i_prb), rec_head.CHSTART };
                    analy.ch_max = rec{ analy.ind_prb(i_prb), rec_head.CHSTOP };
                elseif strcmp(rec{ analy.ind_prb(i_prb), rec_head.SORTDIR }, 'ascending')
                    analy.ch_min = rec{ analy.ind_prb(i_prb), rec_head.CHSTOP };
                    analy.ch_max = rec{ analy.ind_prb(i_prb), rec_head.CHSTART };
                end

                summary.infos{p_ctr} = analy;
                
                infile = [analy.cur_dir '\' analy.monkey '_' num2str(analy.date) '_' ...
                    analy.location '_' num2str(analy.probes) '_June2020_VS' ...
                    '_' summary.d_type{id} '.mat'];
                
                if exist(infile, 'file')
                    load(infile)
                    data_info = eval('info');
                    clear cond info
                    data = T_BLKING_VS(data, data_info);
                    data_cond = T_COND_VS(data, data_info);
                else
                    disp([infile ' \n does not exist']);
                    continue
                end
                
                if strcmp(summary.align_point, 'aot')
                    trial_data = [summary.d_type{id} '_trial'];
                elseif strcmp(summary.align_point, 'aos')
                    trial_data = [summary.d_type{id} '_trial_aos'];
                end
                
                clear both_cond is_cond1 is_cond2 set_size
                
                is_cond1 = data_cond.(summary.conditions{cp, 1});
                is_cond2 = data_cond.(summary.conditions{cp, 2});
                n_cond1 = nansum(is_cond1);
                n_cond2 = nansum(is_cond2);
                both_cond = [find(is_cond1); find(is_cond2)];
                
                dat_n = numel(both_cond);
                if dat_n > 1
                    
                    reli_pr = data.trial_by.block_pres(both_cond);
                    reli_rt = data.trial_by.rt(both_cond);
                    [~,~,rnk_rt] = unique(reli_rt);
                    rnk_rt = rnk_rt./(max(rnk_rt)).*100;

                    reli_loc = abs(data.trial_by.targ_loc(both_cond) -analy.rf_loc +180) -180;                  

                    %switched around vals
                    preave = nanmean(data.(trial_data)(:,1:data_info.evt_pre,both_cond),2);
                    presd = nanstd(data.(trial_data)(:,1:data_info.evt_pre,both_cond),[],2);
                    datlen = size(data.(trial_data)(:,:,both_cond), 2);
                    trllen = size(data.(trial_data)(:,:,both_cond), 3);

                    reli_dat = ((data.(trial_data)(:,:,both_cond) - ...
                        repmat(preave, [1 datlen 1])) ./ ...
                        repmat(presd, [1 datlen 1]));

                    clear preave presd trllen

                    super_chan = analy.L4_bot - 9;
                    deep_chan = analy.L4_bot + 5;
                    if deep_chan > size(reli_dat,1); deep_chan = size(reli_dat,1); end
                    if super_chan < 1; super_chan = 1; end
                    
                    chan_vec = super_chan : deep_chan;
                    
                    for ii = 1:n_cond1

                        reli_dat_pop_targ = cat(1, reli_dat_pop_targ, ...
                            nanmean(reli_dat(chan_vec, data_info.evt_pre-comp_pre:data_info.evt_pre+comp_post, ii), 1), ...
                            reli_dat(chan_vec, data_info.evt_pre-comp_pre:data_info.evt_pre+comp_post, ii));

                        depth_pop_targ = cat(1, depth_pop_targ, (0:numel(chan_vec))');
                        col_pop_targ = cat(1, col_pop_targ, zeros(numel(chan_vec)+1,1)+p_ctr);
                        monk_pop_targ = cat(1, monk_pop_targ, repmat(analy.monkey, numel(chan_vec)+1, 1));
                        rt_targ = cat(1, rt_targ, repmat(reli_rt(ii), numel(chan_vec)+1, 1));
                        pr_targ = cat(1, pr_targ, repmat(reli_pr(ii), numel(chan_vec)+1, 1));
                        rtrnk_targ = cat(1, rtrnk_targ, repmat(rnk_rt(ii), numel(chan_vec)+1, 1));

                    end

                    for ii = 1:n_cond2

                        reli_dat_pop_dist = cat(1, reli_dat_pop_dist, ...
                            nanmean(reli_dat(chan_vec, data_info.evt_pre-comp_pre:data_info.evt_pre+comp_post, n_cond1+ii),1), ...
                            reli_dat(chan_vec, data_info.evt_pre-comp_pre:data_info.evt_pre+comp_post, n_cond1+ii)); %confirmed

                        depth_pop_dist = cat(1, depth_pop_dist, (0:numel(chan_vec))'); %confirmed
                        col_pop_dist = cat(1, col_pop_dist, zeros(numel(chan_vec)+1,1)+p_ctr); %confirmed
                        monk_pop_dist = cat(1, monk_pop_dist, repmat(analy.monkey, numel(chan_vec)+1, 1)); %confirmed
                        rt_dist = cat(1, rt_dist, repmat(reli_rt(ii +n_cond1), numel(chan_vec)+1, 1)); %bad %fixed
                        pr_dist = cat(1, pr_dist, repmat(reli_pr(ii +n_cond1), numel(chan_vec)+1, 1)); %bad %fixed
                        rtrnk_dist = cat(1, rtrnk_dist, repmat(rnk_rt(ii +n_cond1), numel(chan_vec)+1, 1)); %bad %fixed
                        dist_from_targ = cat(1, dist_from_targ, repmat(reli_loc(ii +n_cond1), numel(chan_vec)+1, 1)); %bad %fixed

                    end   
                end
                
                clear reli_rt reli_dat chan_vec super_chan deep_chan n_cond1 n_cond2 dat_n both_cond is_cond1 is_cond2
                clear reli_pr rnk_rt reli_loc
                
            end
        end
        
        save([summary.dir_out summary.d_type{id} '_' ...
            summary.conditions{cp,1} '_' summary.conditions{cp, 2} '_compiled_for_PRAv3.mat'], ...
            '-v7.3', '-nocompression', 'summary', 'reli_dat_pop_targ', ...
            'reli_dat_pop_dist', 'depth_pop_targ', 'depth_pop_dist', ...
            'col_pop_dist', 'col_pop_targ', 'monk_pop_targ', 'monk_pop_dist', ...
            'rt_targ', 'rt_dist', 'pr_targ', 'pr_dist', 'rtrnk_targ', ...
            'rtrnk_dist', 'dist_from_targ');
        
        clear rt_targ rt_dist reli_dat_pop_dist depth_pop_targ reli_dat_pop_targ depth_pop_dist 
        clear col_pop_dist col_pop_targ monk_pop_dist monk_pop_targ dist_from_targ rtrnk_targ rtrnk_dist pr_targ pr_dist
        
    end
end
end