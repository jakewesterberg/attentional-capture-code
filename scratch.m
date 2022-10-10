%%
ucols = unique(col_pop_targ)';
udeps = unique(depth_pop_targ)'; udeps = udeps(2:end);
udat_targ = nan(407, numel(ucols)*numel(udeps), 4);
udat_dist = nan(407, numel(ucols)*numel(udeps), 4);
for i = ucols
    for j = udeps
        for k = 1 : 4

            switch k
                case 1
                    lims = [0 25];
                case 2
                    lims = [25 50];
                case 3
                    lims = [50 75];
                case 4
                    lims = [75 100];
            end

            udat_targ(:,i*j,k) = nanmean(reli_dat_pop_targ( ...
                rtrnk_targ>lims(1) & ...
                rtrnk_targ<=lims(2) &...
                col_pop_targ == i & ...
                depth_pop_targ == j, :)); 

            udat_dist(:,i*j,k) = nanmean(reli_dat_pop_dist( ...
                rtrnk_dist>lims(1) & ...
                rtrnk_dist<=lims(2) &...
                col_pop_dist == i & ...
                depth_pop_dist == j, :)); 

        end
    end
end

%%

ts = -101:305;

% figure
% subplot(1,3,1);  hold on;
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>75 & rtrnk_targ < 101,:)), 'color', [0 1 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>50 & rtrnk_targ < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>25 & rtrnk_targ < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>0 & rtrnk_targ < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)
% 
% subplot(1,3,2);  hold on;
% plot(nanmean(reli_dat_pop_dist(rtrnk_dist>75 & rtrnk_dist < 101,:)), 'color', [0 1 1], 'linewidth', 1.5);
% plot(nanmean(reli_dat_pop_dist(rtrnk_dist>50 & rtrnk_dist < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_dist(rtrnk_dist>25 & rtrnk_dist < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_dist(rtrnk_dist>0 & rtrnk_dist < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)
% 
% subplot(1,3,3);  hold on;
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>75 & rtrnk_targ < 101,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>75 & rtrnk_dist < 101,:)), 'color', [0 1 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>50 & rtrnk_targ < 75,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>50 & rtrnk_dist < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>25 & rtrnk_targ < 50,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>25 & rtrnk_dist < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
% plot(nanmean(reli_dat_pop_targ(rtrnk_targ>0 & rtrnk_targ < 25,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>0 & rtrnk_dist < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)

figure
subplot(1,3,1);  hold on;
plot(nanmean(reli_dat_pop_targ(:,:)), 'color', [0 1 1], 'linewidth', 1.5)
plot(nanmean(reli_dat_pop_targ(dist_from_targ==-20,:)), 'color', [0 .66 1], 'linewidth', 1.5)

subplot(1,3,2);  hold on;
plot(nanmean(reli_dat_pop_dist(:,:)), 'color', [0 1 1], 'linewidth', 1.5);
plot(nanmean(reli_dat_pop_dist(dist_from_targ==-120,:)), 'color', [0 .66 1], 'linewidth', 1.5)

subplot(1,3,3);  hold on;
plot(nanmean(reli_dat_pop_targ(rtrnk_targ>75 & rtrnk_targ < 101,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>75 & rtrnk_dist < 101,:)), 'color', [0 1 1], 'linewidth', 1.5)
plot(nanmean(reli_dat_pop_targ(rtrnk_targ>50 & rtrnk_targ < 75,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>50 & rtrnk_dist < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)


figure; hold on;
for ii = 1 : 4
[ m,l,u ] = confidence_interval(squeeze(udat_targ(:,:,ii)-udat_dist(:,:,ii)));
plot_ci(smooth(l, 25),smooth(u, 25),ts,[0 (ii-1)*.33 1],.15)
plot(ts,smooth(m, 25), 'color', [0 (ii-1)*.33 1], 'linewidth', 1.5)
set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')
end


% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,2)-udat_dist(:,:,2),2),3),10), 'color', [0 .33 1], 'linewidth', 1.5)
% 
% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,3)-udat_dist(:,:,3),2),3),10), 'color', [0 .66 1], 'linewidth', 1.5)
% 
% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,4)-udat_dist(:,:,4),2),3),10), 'color', [0 1 1], 'linewidth', 1.5)
