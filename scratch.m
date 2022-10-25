%%
ucols = unique(col_pop_targ)';
udeps = unique(depth_pop_targ)'; udeps = udeps(2:end);
udat_targ = nan(407, numel(ucols)*numel(udeps), 4);
udat_targ2 = nan(407, numel(ucols)*numel(udeps));
udat_dist = nan(407, numel(ucols)*numel(udeps), 4);
udat_dist_pos = nan(407, numel(ucols)*numel(udeps), 4, 3);
udat_dist_pos2 = nan(407, numel(ucols)*numel(udeps), 4, 3);
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

            udat_targ(:,i*15-15+j,k) = nanmean(reli_dat_pop_targ( ...
                rtrnk_targ>lims(1) & ... 
                rtrnk_targ<=lims(2) &...
                col_pop_targ == i & ...
                depth_pop_targ == j, :)); 

            udat_dist(:,i*15-15+j,k) = nanmean(reli_dat_pop_dist( ...
                rtrnk_dist>lims(1) & ...
                rtrnk_dist<=lims(2) &...
                col_pop_dist == i & ...
                depth_pop_dist == j, :));

            for m = [-60, -120, -180]
                udat_dist_pos(:,i*15-15+j,k,m/-60) = nanmean(reli_dat_pop_dist( ...
                    rtrnk_dist>lims(1) & ...
                    rtrnk_dist<=lims(2) &...
                    col_pop_dist == i & ...
                    dist_from_targ == m & ...
                    depth_pop_dist == j, :));
            end

        end

        udat_targ2(:,i*15-15+j) = nanmean(reli_dat_pop_targ( ...
            col_pop_targ == i & ...
            depth_pop_targ == j, :));

        for m = [-60, -120, -180]
            udat_dist_pos2(:,i*15-15+j,m/-60) = nanmean(reli_dat_pop_dist( ...
                col_pop_dist == i & ...
                dist_from_targ == m & ...
                depth_pop_dist == j, :));
        end
    end
end

%%

ulims = 301:435;
ts = -101:305;

figure
subplot(1,3,1);  hold on;
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>75 & rtrnk_targ < 101,:)), 'color', [0 1 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>50 & rtrnk_targ < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>25 & rtrnk_targ < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>0 & rtrnk_targ < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)
set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')

subplot(1,3,2);  hold on;
plot(ts,nanmean(reli_dat_pop_dist(rtrnk_dist>75 & rtrnk_dist < 101,:)), 'color', [0 1 1], 'linewidth', 1.5);
plot(ts,nanmean(reli_dat_pop_dist(rtrnk_dist>50 & rtrnk_dist < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_dist(rtrnk_dist>25 & rtrnk_dist < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_dist(rtrnk_dist>0 & rtrnk_dist < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)
set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')

subplot(1,3,3);  hold on;
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>75 & rtrnk_targ < 101,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>75 & rtrnk_dist < 101,:)), 'color', [0 1 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>50 & rtrnk_targ < 75,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>50 & rtrnk_dist < 75,:)), 'color', [0 .66 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>25 & rtrnk_targ < 50,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>25 & rtrnk_dist < 50,:)), 'color', [0 .33 1], 'linewidth', 1.5)
plot(ts,nanmean(reli_dat_pop_targ(rtrnk_targ>0 & rtrnk_targ < 25,:)) - nanmean(reli_dat_pop_dist(rtrnk_dist>0 & rtrnk_dist < 25,:)), 'color', [0 0 1], 'linewidth', 1.5)
set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')

% figure
% subplot(1,3,1);  hold on;
% plot(-101:305, nanmean(reli_dat_pop_targ(:,:)), 'color', [0 1 1], 'linewidth', 1.5)
% plot(-101:305, nanmean(reli_dat_pop_dist(dist_from_targ==-60,:)), 'color', [0 .66 1], 'linewidth', 1.5)
% 
% subplot(1,3,2);  hold on;
% plot(-101:305, nanmean(reli_dat_pop_targ(:,:)), 'color', [0 1 1], 'linewidth', 1.5);
% plot(-101:305, nanmean(reli_dat_pop_dist(dist_from_targ==-180,:)), 'color', [0 .66 1], 'linewidth', 1.5)
% 
% subplot(1,3,3);  hold on;
% plot(-101:305, nanmean(reli_dat_pop_targ(:,:)) - nanmean(reli_dat_pop_dist(dist_from_targ==-60,:)), 'color', [0 1 1], 'linewidth', 1.5)
% plot(-101:305, nanmean(reli_dat_pop_targ(:,:)) - nanmean(reli_dat_pop_dist(dist_from_targ==-180,:)), 'color', [0 .66 1], 'linewidth', 1.5)


figure; hold on;
for ii = 1 : 4
    subplot(1,4,ii); hold on;
    [ m,l,u ] = confidence_interval(squeeze(udat_targ(:,ulims,ii)));
    [ m2,l2,u2 ] = confidence_interval(squeeze(udat_dist(:,ulims,ii)));
    plot_ci(smooth(l, 25),smooth(u, 25),ts,[0 (ii-1)*.33 1],.15)
    plot_ci(smooth(l2, 25),smooth(u2, 25),ts,[0 (ii-1)*.33 1],.15)
    plot(ts,smooth(m, 25), 'color', [0 (ii-1)*.33 1], 'linewidth', 1.5)
    plot(ts,smooth(m2, 25), 'color', [0 (ii-1)*.33 1], 'linewidth', 1.5)
    set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')
end

figure; hold on;
for ii = 1 : 4
    [ m,l,u ] = confidence_interval(squeeze(udat_targ(:,ulims,ii)-udat_dist(:,ulims,ii)));
    plot_ci(smooth(l, 25),smooth(u, 25),ts,[0 (ii-1)*.33 1],.15)
    plot(ts,smooth(m, 25), 'color', [0 (ii-1)*.33 1], 'linewidth', 1.5)
    set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')
end

figure; hold on;
for jj = 1 : 3
    subplot(1,3,jj); hold on;
    for ii = 1 : 4
        [ m,l,u ] = confidence_interval(squeeze(udat_targ(:,ulims,ii)-udat_dist_pos(:,ulims,ii,jj)));
        plot_ci(smooth(l, 25),smooth(u, 25),ts,[.2*jj (ii-1)*.33 1],.15)
        plot(ts,smooth(m, 25), 'color', [.2*jj (ii-1)*.33 1], 'linewidth', 1.5)
    end
    set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')
end

figure; hold on;
for ii = 1 : 4
    subplot(1,4,ii); hold on;
    for jj = 1 : 3
        [ m,l,u ] = confidence_interval(squeeze(udat_targ(:,ulims,ii)-udat_dist_pos(:,ulims,ii,jj)));
        plot_ci(smooth(l, 25),smooth(u, 25),ts,[.2*jj (ii-1)*.33 1],.15)
        plot(ts,smooth(m, 25), 'color', [.2*jj (ii-1)*.33 1], 'linewidth', 1.5)
    end
    set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')
end

figure; hold on;
[ m,l,u ] = confidence_interval(squeeze(udat_targ2(:,ulims)));
[ m2,l2,u2 ] = confidence_interval(squeeze(udat_dist_pos2(:,ulims,1)));
[ m3,l3,u3 ] = confidence_interval(squeeze(udat_dist_pos2(:,ulims,2)));
[ m4,l4,u4 ] = confidence_interval(squeeze(udat_dist_pos2(:,ulims,3)));
plot_ci(smooth(l, 25),smooth(u, 25),ts,[0 0 0],.15)
plot_ci(smooth(l2, 25),smooth(u2, 25),ts,[.2*1 (1-1)*.33 1],.15)
plot_ci(smooth(l3, 25),smooth(u3, 25),ts,[.2*2 (2-1)*.33 1],.15)
plot_ci(smooth(l4, 25),smooth(u4, 25),ts,[.2*3 (3-1)*.33 1],.15)
plot(ts,smooth(m, 25), 'color', [0 0 0], 'linewidth', 1.5)
plot(ts,smooth(m2, 25), 'color', [.2*1 (1-1)*.33 1], 'linewidth', 1.5)
plot(ts,smooth(m3, 25), 'color', [.2*2 (2-1)*.33 1], 'linewidth', 1.5)
plot(ts,smooth(m4, 25), 'color', [.2*3 (3-1)*.33 1], 'linewidth', 1.5)
set(gca,'xlim', [-25 200], 'XMinorTick','on','YMinorTick','on','Box','on')

% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,2)-udat_dist(:,:,2),2),3),10), 'color', [0 .33 1], 'linewidth', 1.5)
%
% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,3)-udat_dist(:,:,3),2),3),10), 'color', [0 .66 1], 'linewidth', 1.5)
% 
% plot(ts,smooth(nanmean(nanmean(udat_targ(:,:,4)-udat_dist(:,:,4),2),3),10), 'color', [0 1 1], 'linewidth', 1.5)
