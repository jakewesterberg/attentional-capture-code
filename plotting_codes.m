clear
clc
%close all

% expfit = fittype('b*exp(g*(x-d))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'b', 'g', 'd'});
% 
% pwrfit = fittype('a*(x-d)^b',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a', 'b', 'd'});

% x = ty';
% y = tx';
% 
% pwrfit = @(c)c(1).*(x).^c(2)-y;
% x0 = [1 -1];
% options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',1000000, 'MaxIterations', 100000);
% [z,~,resid,~,~,~,J] = lsqnonlin(pwrfit,x0,[-Inf -Inf], [Inf 0], options);
% ci = nlparci(z,resid,'jacobian', J);
% plot(x,y,'ko',4:250,z(1).*(x).^z(2),'b-'); hold on;
% plot(4:250,ci(1,1).*(x).^ci(2,1),'b-')
% plot(4:250,ci(1,2).*(x).^ci(2,2),'b-'); hold off;



% load('E:\Reliability_Analyses\mua_quarts.mat')
% mua_quarts = out;
% mua_quarts_compiled = compiled;

%figure; subplot1(2,size(out{1}.CS_crit_percent_nrs, 2), 'Gap', [.01 .01])

mult = .75;

% fhc3 = nan(1,mua_quarts_compiled.mpi,4);
% for ii = 5 %: size(mua_quarts{1}.CS_crit_percent_nrs, 2)
% 
%     eval(['p' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     eval(['s' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     eval(['b' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     
%     t_baseline = [];
%     for kk = 1 : numel(mua_quarts)
%         tdat = smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(1:mua_quarts_compiled.mpi,ii,50:100))','movmean',5);
%         %tdat = squeeze(mua_quarts{kk}.CS_crit_percent_nrs(1:mua_quarts_compiled.mpi,ii,50:100));
%         tdat = tdat(:);
%         t_baseline = cat(1, t_baseline, tdat);
%     end
%     sig_thres_up = prctile(t_baseline,100);
%     sig_thres_down = prctile(t_baseline,0);
%     
%     for kk = 1 :numel(mua_quarts)
%         
%         temp_crit_percent_nrs = mua_quarts{kk}.CS_crit_percent_nrs;
%         %temp_crit_percent_nrs = permute(smoothdata(permute(temp_crit_percent_nrs,[3 1 2]), 'movmean', 5),[2 3 1]);
%         
%         temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
%         temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
%         
%         temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
%         temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
%         
%         temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
%         temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
%         
%         temp_exo_nrs_sm = [];
%         temp_endo_nrs_sm = [];
%         
%         n_sm = 49;
%         
%         for iii = abs(mua_quarts_compiled.time_vec(1)) : size(temp_exo_nrs,3)-n_sm
%            
%             temp_exo_nrs_sm(:,:,iii) = sum(temp_exo_nrs(:,:,iii:iii+n_sm),3);
%             temp_endo_nrs_sm(:,:,iii) = sum(temp_endo_nrs(:,:,iii:iii+n_sm),3);
%             
%         end
% 
%         for jj = 1 : mua_quarts_compiled.mpi
%             
%             if kk == 4
%                 mrk1 = [0,1,1];
%                 mrk2 = [1*mult,0,0];
%             elseif kk == 3
%                 mrk1 = [0,.66,1];
%                 mrk2 = [1*mult,.33*mult,0];
%             elseif kk == 2
%                 mrk1 = [0,.33,1];
%                 mrk2 = [1*mult,.66*mult,0];
%             else
%                 mrk1 = [0,0,1];
%                 mrk2 = [1*mult,1*mult,0];
%             end
%             
%             set(0,'CurrentFigure',eval(['p' num2str(ii)])); hold on;
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_quarts_compiled.mpi)).*.66);
% 
%             if jj == 1 && kk == 1
%                 
%                 set(gca, 'linewidth', 2, 'ylim', [1 100], 'xlim', [-100 225], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
%                 patch([-100 0 0 -100], [0 0 33.332 33.332], ...
%                     'k', 'edgecolor', 'none',  'facealpha', .1)
%                 
%                 vt = plot([0 0], [0 33.322+8.333]);
%                 set(vt, 'color', [0 0 0 .75], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([58 78 78 58], ...
%                     [0 0 100 100], ...
%                     [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 vt = vline([58 78]);
%                 set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([mua_quarts_compiled.time_vec(1) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(1)], ...
%                     [sig_thres_down sig_thres_down sig_thres_up sig_thres_up], ...
%                     'red', 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 ht = hline([sig_thres_up sig_thres_down]);
%                 set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%             end
%             
%             if jj == 1
%                 
% %                 patch([round(nanmean(mua_quarts{kk}.rt_med(mua_quarts_compiled.mpi,ii,:)))-10 round(nanmean(mua_quarts{kk}.rt_med(mua_quarts_compiled.mpi,ii,:))) ...
% %                     round(nanmean(mua_quarts{kk}.rt_med(mua_quarts_compiled.mpi,ii,:))) round(nanmean(mua_quarts{kk}.rt_med(mua_quarts_compiled.mpi,ii,:)))- 10],...
% %                     [0 0 100 100], ...
% %                     mrk1, 'edgecolor', 'none',  'facealpha', .25)
%                 
%             end
%             
%             if  jj == mua_quarts_compiled.mpi
%                 
% %                 plot([round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:))) round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:)))], ...
% %                     [0 100], 'color', mrk1, 'linewidth', 2)
%                 
%             end
%             
%             endingp = find(mua_quarts_compiled.time_vec==round(((nanmean(mua_quarts{kk}.rt_med(jj,ii,:))))-10));
%             if isempty(endingp); endingp = 301; end
%             
%             %find(mua_quarts{kk}.pcttotl(1,ii,:) < 99,1)-10)
%             plot(mua_quarts_compiled.time_vec(1:endingp), ...
%                 ...smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:301)),'movmean',5),...
%                 squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)), ...
%                 'Color', [tmrk .5], 'linewidth', 2);
%             
%             
%             if kk == numel(mua_quarts) && jj == mua_quarts_compiled.mpi
%                 ht = hline(16.666);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [-50 round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:)))])
%             end
%             
%             %vt = vline(round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:))));
%             %set(vt,'linewidth',1,'linestyle', '-', 'color',[mrk1 1/250])
%             
%             fhc = find(squeeze(temp_hit_crit_nrs(jj,ii,:)),1);
%             fhc2 = find(squeeze(temp_exo_nrs_sm(jj,ii,:))==(n_sm+1),1);
%             
%             %subplot1(ii+size(mua_quarts{kk}.CS_crit_percent_nrs, 2)); hold on;
%             set(0,'CurrentFigure',eval(['s' num2str(ii)])); hold on;
%             
%             if ~isempty(fhc)
%                 if fhc+mua_quarts_compiled.time_vec(1) < round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:)))
%                     scatter(nanmean(mua_quarts{kk}.lat_5(jj,ii,:)), mua_quarts_compiled.counts(jj), 'markerfacecolor', ...
%                         tmrk, 'markeredgecolor', [0 0 0]);
%                     %scatter(fhc+mua_quarts_compiled.time_vec(1), mua_quarts_compiled.counts(jj), 'markerfacecolor', ...
%                     %    tmrk, 'markeredgecolor', [0 0 0]); clear fhc
%                     scatter(fhc2+mua_quarts_compiled.time_vec(1), mua_quarts_compiled.counts(jj), 'markerfacecolor', ...
%                         tmrk, 'markeredgecolor', [0 0 0]);
%                     set(gca, 'linewidth', 2, 'ylim', [1 mua_quarts_compiled.mpi], ...
%                         'xlim', [-50 250], 'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                     box off
%                 end
%             end
%             if jj == 1
%                 vt = vline(0);
%                 set(vt, 'color', 'k', 'linewidth', 2, 'linestyle', ':')
%             end
%             
%             if  jj == mua_quarts_compiled.mpi && kk == 1
%                 set(gca, 'xlim', [-50 300])
%             end
%             
%             
%             
%             set(0,'CurrentFigure',eval(['b' num2str(ii)])); hold on;
%             
%             if jj == 1 && kk == 1
%                 
%                 patch([0 20 20 0], ...
%                     [0 0 250 250], ...
%                     [1 .4 0], 'edgecolor', 'k',  'facealpha', .3)
%                 
%                 patch([0 -3 -3 0], ...
%                     [0 0 250 250], ...
%                     [1 .4 0], 'edgecolor', 'k',  'facealpha', .15)
%      
%             end
% 
%             if ~isempty(fhc2)
%                 fhc3(ii,jj,kk) = fhc2+mua_quarts_compiled.time_vec(1)-nanmedian(mua_quarts{kk}.lat_5(jj,ii,:));
%                 %if fhc+mua_quarts_compiled.time_vec(1) < round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:)))
%                     
%                     %plot([(fhc2+mua_quarts_compiled.time_vec(1)-nanmean(mua_quarts{kk}.lat_5(jj,ii,:))) 300], [jj-kk/5 jj-kk/5], ...
%                     %    'color', mrk1, 'linewidth', 1)
%                     
%                      scatter(fhc2+mua_quarts_compiled.time_vec(1)-nanmedian(mua_quarts{kk}.lat_5(jj,ii,:)), jj, ...
%                          'markerfacecolor', tmrk, 'markeredgecolor', mrk1, 'linewidth', 1)
%                      
%                      set(gca, 'linewidth', 2, 'ylim', [0 250], ...
%                          'xlim', [-10 175], 'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                     box off
%                     
%                 %end
%             end
%             if jj == 1
%                 vt = vline(0);
%                 set(vt, 'color', 'k', 'linewidth', 2, 'linestyle', ':')
%             end
%             
%             if  jj == mua_quarts_compiled.mpi && kk == 1
%                 set(gca, 'xlim', [-50 round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:)))])
%             end
%             
%             if jj == mua_quarts_compiled.mpi
%                 
%                  tx = squeeze(fhc3(ii,1:end,kk)+mua_quarts_compiled.time_vec(1))' - squeeze(nanmean(mua_quarts{kk}.lat_5(1:mua_quarts_compiled.mpi,ii,:)+2,3));
%                  ty = (1:mua_quarts_compiled.mpi)';
%                  
%                  ty=ty(~isnan(tx));
%                  tx=tx(~isnan(tx));
% %                 
% %                 tfit = fit(tx, ty, pwrfit, ...
% %                      'Lower', [-Inf -1 -Inf], ...
% %                      'Upper', [Inf -.00001 Inf]);
% %                 
% % %                 tfit = fit(tx, ty, expfit, ...
% % %                     'Lower', [0 -Inf], ...
% % %                     'Upper', [Inf 0]);
% %                 
% %                 plot(tfit, tx, ty, 'predfunc', .95); legend off
%                 
%             end
%         end
%     end
% end
% 
% close all
% figure('Renderer', 'Painters');
% subplot1(15,1);
% fhc3 = nan(1,mua_quarts_compiled.mpi,4);
% for ii = 6 : size(mua_quarts{1}.CS_crit_percent_nrs, 2)
% 
%     subplot1(ii-5)
%     
%     t_baseline = [];
%     for kk = 1 : numel(mua_quarts)
%         tdat = smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(1:mua_quarts_compiled.mpi,ii,50:100))','movmean',10);
%         %tdat = squeeze(mua_quarts{kk}.CS_crit_percent_nrs(1:mua_quarts_compiled.mpi,ii,50:100));
%         tdat = tdat(:);
%         t_baseline = cat(1, t_baseline, tdat);
%     end
%     sig_thres_up = prctile(t_baseline,99);
%     sig_thres_down = prctile(t_baseline,1);
%     
%     for kk = [1 4]%: numel(mua_quarts)
%         
%         temp_crit_percent_nrs = mua_quarts{kk}.CS_crit_percent_nrs;
%         %temp_crit_percent_nrs = permute(smoothdata(permute(temp_crit_percent_nrs,[3 1 2]), 'movmean', 5),[2 3 1]);
%         
%         temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
%         temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
%         
%         temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
%         temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
%         
%         temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
%         temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
%         
%         n_sm = 49;
%        
%         for jj =  mua_quarts_compiled.mpi
%             
%             if kk == 4
%                 mrk1 = [0,1,1];
%                 mrk2 = [1*mult,0,0];
%             elseif kk == 3
%                 mrk1 = [0,.66,1];
%                 mrk2 = [1*mult,.33*mult,0];
%             elseif kk == 2
%                 mrk1 = [0,.33,1];
%                 mrk2 = [1*mult,.66*mult,0];
%             else
%                 mrk1 = [0,0,1];
%                 mrk2 = [1*mult,1*mult,0];
%             end
%             
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_quarts_compiled.mpi)).*.66);
% 
%             if jj == 250 && kk == 1
%                 
%                 set(gca, 'linewidth', 2, 'ylim', [0 100], 'xlim', [-100 200], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
%                 if ii == 20
%                     patch([58 78 78 58], ...
%                         [0 0 100 100], ...
%                         [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
%                     
%                     
%                     vt = vline([58 78]);
%                     set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 end
%                 
%                 
%                 
%                 patch([mua_quarts_compiled.time_vec(1) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(1)], ...
%                     [sig_thres_down sig_thres_down sig_thres_up sig_thres_up], ...
%                     'red', 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 ht = hline([sig_thres_up sig_thres_down]);
%                 set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%             end
%             
%             endingp = 251;
%             
%             plot(mua_quarts_compiled.time_vec(1:endingp), ...
%                 smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5),...
%                 ...squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)), ...
%                 'Color', [tmrk .5], 'linewidth', 2);
%             
%             fill([mua_quarts_compiled.time_vec(1:endingp) fliplr(mua_quarts_compiled.time_vec(1:endingp))], ...
%                 [zeros(1,251)+16.666, ...
%                 fliplr(smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5)')],...
%                 ...squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)), ...
%                 [tmrk], 'edgecolor', 'none', 'facealpha', 1);
%             
%             
%             if kk == numel(mua_quarts) && jj == mua_quarts_compiled.mpi
%                 ht = hline(16.666);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [0 200])
%             end
%             
%             %vt = vline(round(nanmean(mua_quarts{kk}.rt_med(jj,ii,:))));
%             %set(vt,'linewidth',1,'linestyle', '-', 'color',[mrk1 1/250])
%             
%         end
%     end
%     
%     if (ii - 5)  < 15
%         
%         axis off
%         
%     end
%             
% end






% load('E:\Reliability_Analyses\mua_quarts_aos.mat')
% mua_quarts_aos = out;
% mua_quarts_aos_compiled = compiled;
% close all
% figure('Renderer', 'Painters'); hold on;
% fhc3 = nan(1,mua_quarts_aos_compiled.mpi,4);
% for ii = 5 %: size(mua_quarts{1}.CS_crit_percent_nrs, 2)
% 
%     for kk = 1 : numel(mua_quarts_aos)
%         
%         n_sm = 9;
%        
%         for jj =  1:mua_quarts_aos_compiled.mpi
%             
%             if kk == 4
%                 mrk1 = [0,1,1];
%                 mrk2 = [1*mult,0,0];
%             elseif kk == 3
%                 mrk1 = [0,.66,1];
%                 mrk2 = [1*mult,.33*mult,0];
%             elseif kk == 2
%                 mrk1 = [0,.33,1];
%                 mrk2 = [1*mult,.66*mult,0];
%             else
%                 mrk1 = [0,0,1];
%                 mrk2 = [1*mult,1*mult,0];
%             end
%             
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_quarts_aos_compiled.mpi)).*.66);
% 
%             if jj == 250 && kk == 1
%                 
%                 set(gca, 'linewidth', 2, 'ylim', [0 100], 'xlim', [-100 10], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
%             end
%             
%             plot(mua_quarts_aos_compiled.time_vec(50:161), ...
%                 ...smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5),...
%                 squeeze(mua_quarts_aos{kk}.CS_crit_percent_nrs(jj,ii,50:161)), ...
%                 'Color', [tmrk .5], 'linewidth', 2);
%             
%             if kk == numel(mua_quarts_aos) && jj == mua_quarts_aos_compiled.mpi
%                 ht = hline(16.666);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [-75 10])
%             end
%             
%         end
%     end
% end
% 
% 




% clear varfhc hvar pvar
% for j = 1 : 4
%     for i = 1 : 250
%         
%         if i > 240
%             varfhc(i,j) = var(fhc3(241:250,j));
%         else
%             varfhc(i,j) = var(fhc3(i:i+9,j));
%         end
%         
%         %[hvar(i,j), pvar(i,j)] = vartest2(fhc3(i:i+14,j), fhc3(i+1:i+15,j));
%         
%     end
% end
% 
% load('E:\Reliability_Analyses\mua_priming.mat')
% mua_priming = out;
% mua_priming_compiled = compiled;
% 
% clear out compiled
% 
% close all
% mua_priming_compiled.mpi = 250;
% 
% mult = .75;
% for ii = [5] %size(mua_priming{1}.CS_crit_percent_nrs, 2)
% 
%     eval(['p' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     eval(['s' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     
%     t_baseline = [];
%     for kk = 1 : numel(mua_priming)
%         tdat = squeeze(mua_priming{kk}.CS_crit_percent_nrs(1:mua_priming_compiled.mpi,ii,1:100));
%         tdat = tdat(:);
%         t_baseline = cat(1, t_baseline, tdat);
%     end
%     sig_thres_up = prctile(t_baseline,99);
%     sig_thres_down = prctile(t_baseline,1);
%     
%     for kk = 2 %numel(mua_priming):-1:1
%         
%         temp_crit_percent_nrs = mua_priming{kk}.CS_crit_percent_nrs;
%         
%         temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
%         temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
%         
%         temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
%         temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
%         
%         temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
%         temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
%         
%         temp_exo_nrs_sm = [];
%         temp_endo_nrs_sm = [];
%         
%         n_sm = 4;
%         
%         for iii = abs(mua_priming_compiled.time_vec(1)) : size(temp_exo_nrs,3)-n_sm
%            
%             temp_exo_nrs_sm(:,:,iii) = sum(temp_exo_nrs(:,:,iii:iii+n_sm),3);
%             temp_endo_nrs_sm(:,:,iii) = sum(temp_endo_nrs(:,:,iii:iii+n_sm),3);
%             
%         end
% 
%         for jj = 1 : mua_priming_compiled.mpi
%             
%             if kk == 1
%                 mrk1 = [1,0,1];
%                 mrk2 = [0,1*mult,0];
%             else
%                 mrk1 = [0,0,0];
%                 mrk2 = [1*mult,1*mult,1*mult];
%             end
%             
%             set(0,'CurrentFigure',eval(['p' num2str(ii)])); hold on;
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_priming_compiled.mpi)).*.66);
% 
%             if jj == 1 && kk == numel(mua_priming)
%                 
%                 set(gca, 'linewidth', 2, 'ylim', [1 100], 'xlim', [-100 300], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
% %                 patch([-100 0 0 -100], [0 0 33.332+8.333 33.332+8.333], ...
% %                     'k', 'edgecolor', 'none',  'facealpha', .1)
%                 
% %                 vt = plot([0 0], [0 33.322+16.666]);
% %                 set(vt, 'color', [0 0 0 .75], 'linewidth', 2, 'linestyle', '-')
%                 
% %                 patch([53 78 78 53], ...
% %                     [0 0 100 100], ...
% %                     [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
% %                 
% %                 vt = vline([53 78]);
% %                 set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([mua_priming_compiled.time_vec(1) mua_priming_compiled.time_vec(end) mua_priming_compiled.time_vec(end) mua_priming_compiled.time_vec(1)], ...
%                     [sig_thres_down sig_thres_down sig_thres_up sig_thres_up], ...
%                     'red', 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 ht = hline([sig_thres_up sig_thres_down]);
%                 set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%             end
%             
%             if jj == 1
%                 
% %                 patch([round(nanmean(mua_priming{kk}.rt_med(mua_priming_compiled.mpi,ii,:)))-10 round(nanmean(mua_priming{kk}.rt_med(mua_priming_compiled.mpi,ii,:))) ...
% %                     round(nanmean(mua_priming{kk}.rt_med(mua_priming_compiled.mpi,ii,:))) round(nanmean(mua_priming{kk}.rt_med(mua_priming_compiled.mpi,ii,:)))- 10],...
% %                     [0 0 100 100], ...
% %                     mrk1, 'edgecolor', 'none',  'facealpha', .25)
%                 
%             end
%             
%             if  jj == mua_priming_compiled.mpi
%                 
% %                 plot([round(nanmean(mua_priming{kk}.rt_med(jj,ii,:))) round(nanmean(mua_priming{kk}.rt_med(jj,ii,:)))], ...
% %                     [0 100], 'color', mrk1, 'linewidth', 2)
%                 
%             end
%             
%             endingp = find(mua_priming_compiled.time_vec==round((nanmean(mua_priming{kk}.rt_med(jj,ii,:)))));
%             if isempty(endingp); endingp = 201; end
%             plot(mua_priming_compiled.time_vec(1:endingp), ...
%                 squeeze(mua_priming{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)), ...
%                 'Color', tmrk, 'linewidth', 2);
%             
%             if kk == 1 && jj == mua_priming_compiled.mpi
%                 ht = hline(16.666);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [-50 round(nanmean(mua_priming{kk}.rt_med(jj,ii,:)))])
%             end
%             
%             %vt = vline(round(nanmean(mua_priming{kk}.rt_med(jj,ii,:))));
%             %set(vt,'linewidth',1,'linestyle', '-', 'color',[mrk1 1/250])
%             
%             fhc = find(squeeze(temp_hit_crit_nrs(jj,ii,:)),1);
%             fhc2 = find(squeeze(temp_exo_nrs_sm(jj,ii,:))==(n_sm+1),1);
%             
% %             subplot1(ii+size(mua_priming{kk}.CS_crit_percent_nrs, 2)); hold on;
%             set(0,'CurrentFigure',eval(['s' num2str(ii)])); hold on;
%             
%             if ~isempty(fhc)
%                 if fhc+mua_priming_compiled.time_vec(1) < round(nanmean(mua_priming{kk}.rt_med(jj,ii,:)))
%                     scatter(nanmean(mua_priming{kk}.lat_5(jj,ii,:)), mua_priming_compiled.counts(jj), 'markerfacecolor', ...
%                         tmrk, 'markeredgecolor', [0 0 0]);
%                     scatter(fhc+mua_priming_compiled.time_vec(1), mua_priming_compiled.counts(jj), 'markerfacecolor', ...
%                         tmrk, 'markeredgecolor', [0 0 0]); clear fhc
%                     scatter(fhc2+mua_priming_compiled.time_vec(1), mua_priming_compiled.counts(jj), 'markerfacecolor', ...
%                         tmrk, 'markeredgecolor', [0 0 0]); clear fhc2
%                     set(gca, 'linewidth', 2, 'ylim', [1 mua_priming_compiled.mpi], ...
%                         'xlim', [-50 250], 'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                     box off
%                 end
%             end
%             if jj == 1
%                 vt = vline(0);
%                 set(vt, 'color', 'k', 'linewidth', 2, 'linestyle', ':')
%             end
%             
%             if  jj == mua_priming_compiled.mpi && kk == 1
%                 set(gca, 'xlim', [-50 round(nanmean(mua_priming{kk}.rt_med(jj,ii,:)))])
%             end
%         end
%     end
% end
% 

% load('E:\Reliability_Analyses\mua_err2.mat')
% mua_err = out;
% mua_err_compiled = compiled;
% 
% clear out compiled
% 
% mult = .75;
% 
% %close all
% mua_err_compiled.mpi = 250;
% 
% fhc3 = nan(1,mua_err_compiled.mpi,4);
% for ii = [5]%: size(mua_quarts{1}.CS_crit_percent_nrs, 2)
% 
%     eval(['p' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     %eval(['s' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     %eval(['b' num2str(ii) '=figure(''Renderer'',''Painters'');'])
%     
%     t_baseline = [];
%     for kk = 1 : numel(mua_err)
%         tdat = smoothdata(squeeze(mua_err{kk}.CS_crit_percent_nrs(1:mua_err_compiled.mpi,ii,50:100))','movmean',5);
%         tdat = squeeze(mua_err{kk}.CS_crit_percent_nrs(1:mua_err_compiled.mpi,ii,50:100));
%         tdat = tdat(:);
%         t_baseline = cat(1, t_baseline, tdat);
%     end
%     sig_thres_up = prctile(t_baseline,99.9);
%     sig_thres_down = prctile(t_baseline,.1);
%     
%     for kk = 1 :numel(mua_err)
%         
%         temp_crit_percent_nrs = mua_err{kk}.CS_crit_percent_nrs;
%         temp_crit_percent_nrs = permute(smoothdata(permute(temp_crit_percent_nrs,[3 1 2]), 'movmean', 5),[2 3 1]);
%         
%         temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
%         temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
%         
%         temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
%         temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
%         
%         temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
%         temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
%         
%         temp_exo_nrs_sm = [];
%         temp_endo_nrs_sm = [];
%         
%         n_sm = 9;
%         
%         for iii = abs(mua_err_compiled.time_vec(1)) : size(temp_exo_nrs,3)-n_sm
%            
%             temp_exo_nrs_sm(:,:,iii) = sum(temp_exo_nrs(:,:,iii:iii+n_sm),3);
%             temp_endo_nrs_sm(:,:,iii) = sum(temp_endo_nrs(:,:,iii:iii+n_sm),3);
%             
%         end
% 
%         for jj = 1 : mua_err_compiled.mpi
%             
%             if kk == 4
%                 mrk1 = [0,1,1];
%                 mrk2 = [1*mult,0,0];
%             elseif kk == 3
%                 mrk1 = [0,.66,1];
%                 mrk2 = [1*mult,.33*mult,0];
%             elseif kk == 2
%                 mrk1 = [0,.33,1];
%                 mrk2 = [1*mult,.66*mult,0];
%             else
%                 mrk1 = [1,0,.75];
%                 mrk2 = [0,1*mult,.25*mult];
%             end
%             
%             set(0,'CurrentFigure',eval(['p' num2str(ii)])); hold on;
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_err_compiled.mpi)).*.66);
% 
%             if jj == 1 && kk == 1
%                 
%                 set(gca, 'linewidth', 2, 'ylim', [-15 15], 'xlim', [-50 250], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
%                 %patch([-100 0 0 -100], [0 0 33.332+8.333 33.332+8.333], ...
%                 %    'k', 'edgecolor', 'none',  'facealpha', .1)
%                 
%                 vt = plot([0 0], [-100 100]);
%                 set(vt, 'color', [0 0 0 .75], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([58 78 78 58], ...
%                     [-100 -100 100 100], ...
%                     [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 vt = vline([58 78]);
%                 set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([mua_err_compiled.time_vec(1) mua_err_compiled.time_vec(end) mua_err_compiled.time_vec(end) mua_err_compiled.time_vec(1)], ...
%                     [sig_thres_down-16.666 sig_thres_down-16.666 sig_thres_up-16.666 sig_thres_up-16.666], ...
%                     'red', 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 ht = hline([sig_thres_up-16.666 sig_thres_down-16.666]);
%                 set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%             end
%             
%             if jj == 1
%                 
% %                 patch([round(nanmean(mua_err{kk}.rt_med(mua_err_compiled.mpi,ii,:)))-10 round(nanmean(mua_err{kk}.rt_med(mua_err_compiled.mpi,ii,:))) ...
% %                     round(nanmean(mua_err{kk}.rt_med(mua_err_compiled.mpi,ii,:))) round(nanmean(mua_err{kk}.rt_med(mua_err_compiled.mpi,ii,:)))- 10],...
% %                     [0 0 100 100], ...
% %                     mrk1, 'edgecolor', 'none',  'facealpha', .25)
%                 
%             end
%             
%             if  jj == mua_err_compiled.mpi
%                 
% %                 plot([round(nanmean(mua_err{kk}.rt_med(jj,ii,:))) round(nanmean(mua_err{kk}.rt_med(jj,ii,:)))], ...
% %                     [0 100], 'color', mrk1, 'linewidth', 2)
%                 
%             end
%             
%             endingp = find(mua_err_compiled.time_vec==round((nanmean(mua_err{kk}.rt_med(jj,ii,:)))));
%             if isempty(endingp); endingp = 301; end
%             plot(mua_err_compiled.time_vec(1:endingp), ...
%                 smoothdata(squeeze(mua_err{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5)-16.666,...
%                 'Color', [tmrk .5], 'linewidth', 2);
%             
%             
%             if kk == numel(mua_err) && jj == mua_err_compiled.mpi
%                 ht = hline(0);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [-50 250])
%             end
%             
%             %vt = vline(round(nanmean(mua_err{kk}.rt_med(jj,ii,:))));
%             %set(vt,'linewidth',1,'linestyle', '-', 'color',[mrk1 1/250])
%             
%             fhc = find(squeeze(temp_hit_crit_nrs(jj,ii,:)),1);
%             fhc2 = find(squeeze(temp_endo_nrs_sm(jj,ii,:))==(n_sm+1),1);
%             %fhc3(ii,1:end,kk)+mua_quarts_compiled.time_vec(1))' - squeeze(nanmean(mua_quarts{kk}.lat_5(1:mua_quarts_compiled.mpi,ii,:),3)
%             
%             %subplot1(ii+size(mua_err{kk}.CS_crit_percent_nrs, 2)); hold on;
% %             set(0,'CurrentFigure',eval(['s' num2str(ii)])); hold on;
% %             
% %             if ~isempty(fhc)
% %                 if fhc+mua_err_compiled.time_vec(1) < round(nanmean(mua_err{kk}.rt_med(jj,ii,:)))
% %                     scatter(nanmean(mua_err{kk}.lat_5(jj,ii,:)), mua_err_compiled.counts(jj), 'markerfacecolor', ...
% %                         tmrk, 'markeredgecolor', [0 0 0]);
% %                     %scatter(fhc+mua_err_compiled.time_vec(1), mua_err_compiled.counts(jj), 'markerfacecolor', ...
% %                     %    tmrk, 'markeredgecolor', [0 0 0]); clear fhc
% %                     scatter(fhc2+mua_err_compiled.time_vec(1), mua_err_compiled.counts(jj), 'markerfacecolor', ...
% %                         tmrk, 'markeredgecolor', [0 0 0]);
% %                     set(gca, 'linewidth', 2, 'ylim', [1 mua_err_compiled.mpi], ...
% %                         'xlim', [-50 300], 'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
% %                     box off
% %                 end
% %             end
% %             if jj == 1
% %                 vt = vline(0);
% %                 set(vt, 'color', 'k', 'linewidth', 2, 'linestyle', ':')
% %             end
% %             
% %             if  jj == mua_err_compiled.mpi && kk == 1
% %                 set(gca, 'xlim', [-50 225])
% %             end
% %             
% %             
% %             
% %             set(0,'CurrentFigure',eval(['b' num2str(ii)])); hold on;
% %             
% %             if jj == 1 && kk == 1
% %                 
% %                 patch([0 100 100 0], ...
% %                     [0 0 20 20], ...
% %                     [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
% %      
% %             end
% 
%             if ~isempty(fhc2)
%                 fhc3(ii,jj,kk) = fhc2;
%             end
% 
%         end
%     end
% end
% 
% 
% 
% 
% 
% 
% load('E:\Reliability_Analyses\mua_err2.mat')
% mua_err = out;
% mua_err_compiled = compiled;
% 
% clear out compiled
% 
% mult = .75;
% 
% %close all
% mua_err_compiled.mpi = 250;
% 
% fhc3 = nan(1,mua_err_compiled.mpi,4);
% 
% figure
% subplot1(3,1);
% 
% 
% for ii = [2:4]%: size(mua_quarts{1}.CS_crit_percent_nrs, 2)
%     
%     subplot1(ii-1)
%     t_baseline = [];
%     for kk = 1 : numel(mua_err)
%         tdat = smoothdata(squeeze(mua_err{kk}.CS_crit_percent_nrs(1:mua_err_compiled.mpi,ii,1:50))','movmean',5);
%         tdat = squeeze(mua_err{kk}.CS_crit_percent_nrs(1:mua_err_compiled.mpi,ii,1:50));
%         tdat = tdat(:);
%         t_baseline = cat(1, t_baseline, tdat);
%     end
%     sig_thres_up = prctile(t_baseline,99.9);
%     sig_thres_down = prctile(t_baseline,.1);
%     
%     for kk = 1 :numel(mua_err)
%         
%         temp_crit_percent_nrs = mua_err{kk}.CS_crit_percent_nrs;
%         temp_crit_percent_nrs = permute(smoothdata(permute(temp_crit_percent_nrs,[3 1 2]), 'movmean', 5),[2 3 1]);
%         
%         temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
%         temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
%         
%         temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
%         temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
%         
%         temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
%         temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
%         
%         temp_exo_nrs_sm = [];
%         temp_endo_nrs_sm = [];
%         
%         n_sm = 9;
%         
%         for iii = abs(mua_err_compiled.time_vec(1)) : size(temp_exo_nrs,3)-n_sm
%            
%             temp_exo_nrs_sm(:,:,iii) = sum(temp_exo_nrs(:,:,iii:iii+n_sm),3);
%             temp_endo_nrs_sm(:,:,iii) = sum(temp_endo_nrs(:,:,iii:iii+n_sm),3);
%             
%         end
% 
%         for jj = 1 : mua_err_compiled.mpi
%             
%             if kk == 4
%                 mrk1 = [0,1,1];
%                 mrk2 = [1*mult,0,0];
%             elseif kk == 3
%                 mrk1 = [0,.66,1];
%                 mrk2 = [1*mult,.33*mult,0];
%             elseif kk == 2
%                 mrk1 = [0,.33,1];
%                 mrk2 = [1*mult,.66*mult,0];
%             else
%                 mrk1 = [1,0,.75];
%                 mrk2 = [0,1*mult,.25*mult];
%             end
% 
%             tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_err_compiled.mpi)).*.66);
% 
%             if jj == 1 && kk == 1
%                 
%               set(gca, 'linewidth', 2, 'ylim', [-25 25], 'xlim', [-50 250], ...
%                     'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
%                 box off
%                 
%                 %patch([-100 0 0 -100], [0 0 33.332+8.333 33.332+8.333], ...
%                 %    'k', 'edgecolor', 'none',  'facealpha', .1)
%                 
%                 vt = plot([0 0], [-100 100]);
%                 set(vt, 'color', [0 0 0 .75], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([58 78 78 58], ...
%                     [-100 -100 100 100], ...
%                     [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 vt = vline([58 78]);
%                 set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%                 patch([mua_err_compiled.time_vec(1) mua_err_compiled.time_vec(end) mua_err_compiled.time_vec(end) mua_err_compiled.time_vec(1)], ...
%                     [sig_thres_down-16.666 sig_thres_down-16.666 sig_thres_up-16.666 sig_thres_up-16.666], ...
%                     'red', 'edgecolor', 'none',  'facealpha', .25)
%                 
%                 ht = hline([sig_thres_up-16.666 sig_thres_down-16.666]);
%                 set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
%                 
%                 
%             end
% 
%             
%             endingp = 301;
%             plot(mua_err_compiled.time_vec(1:endingp), ...
%                 smoothdata(squeeze(mua_err{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5)-16.666,...
%                 'Color', [tmrk .5], 'linewidth', 2);
%             
%             
%             if kk == numel(mua_err) && jj == mua_err_compiled.mpi
%                 ht = hline(0);
%                 set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
%                 set(gca, 'xlim', [-50 250])
%             end
%             
%             %vt = vline(round(nanmean(mua_err{kk}.rt_med(jj,ii,:))));
%             %set(vt,'linewidth',1,'linestyle', '-', 'color',[mrk1 1/250])
%             
%             fhc = find(squeeze(temp_hit_crit_nrs(jj,ii,:)),1);
%             fhc2 = find(squeeze(temp_endo_nrs_sm(jj,ii,:))==(n_sm+1),1);
% 
%             if ~isempty(fhc2)
%                 fhc3(ii,jj,kk) = fhc2;
%             end
% 
%         end
%     end
% end













tident='hate_quarts';
load(['E:\Reliability_Analyses\' tident '.mat'])
mua_quarts = out;
mua_quarts_compiled = compiled;
mua_quarts_compiled.mpi = 250;

for ii = 5 %: size(mua_quarts{1}.CS_crit_percent_nrs, 2)

    eval(['p' num2str(ii) '=figure(''Renderer'',''Painters'');'])
    
    t_baseline = [];
    for kk = [1 4] %1 : numel(mua_quarts)
        tdat = smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(1:mua_quarts_compiled.mpi,ii,1:100))','movmean',5);
        tdat = tdat(:);
        t_baseline = cat(1, t_baseline, tdat);
    end
    sig_thres_up = prctile(t_baseline,99.5);
    sig_thres_down = prctile(t_baseline,.5);
    
    for kk = [1 4] %:numel(mua_quarts)
        
        temp_crit_percent_nrs = mua_quarts{kk}.CS_crit_percent_nrs;
        temp_crit_percent_nrs = permute(smoothdata(permute(temp_crit_percent_nrs,[3 1 2]), 'movmean', 5),[2 3 1]);
        
        temp_hit_crit_nrs = temp_crit_percent_nrs == 100;
        temp_hit_crit_nrs(isnan(temp_hit_crit_nrs)) = 0;
        
        temp_exo_nrs = temp_crit_percent_nrs > sig_thres_up;
        temp_exo_nrs(isnan(temp_exo_nrs)) = 0;
        
        temp_endo_nrs = temp_crit_percent_nrs < sig_thres_down;
        temp_endo_nrs(isnan(temp_endo_nrs)) = 0;
        
        temp_exo_nrs_sm = [];
        temp_endo_nrs_sm = [];
        
        n_sm = 49;
        
        for iii = abs(mua_quarts_compiled.time_vec(1)) : size(temp_exo_nrs,3)-n_sm
           
            temp_exo_nrs_sm(:,:,iii) = sum(temp_exo_nrs(:,:,iii:iii+n_sm),3);
            temp_endo_nrs_sm(:,:,iii) = sum(temp_endo_nrs(:,:,iii:iii+n_sm),3);
            
        end

        for jj = 1 : mua_quarts_compiled.mpi
            
            if kk == 4
                mrk1 = [0,1,1];
                mrk2 = [1*mult,0,0];
                endingp = 347;
            elseif kk == 3
                mrk1 = [0,.66,1];
                mrk2 = [1*mult,.33*mult,0];
                endingp = 305;
            elseif kk == 2
                mrk1 = [0,.33,1];
                mrk2 = [1*mult,.66*mult,0];
                endingp = 275;
            else
                mrk1 = [0,0,1];
                mrk2 = [1*mult,1*mult,0];
                endingp = 251;
            end
            
            set(0,'CurrentFigure',eval(['p' num2str(ii)])); hold on;
            tmrk = [1 1 1] - (mrk2.*.33) - ((mrk2.*(jj / mua_quarts_compiled.mpi)).*.66);

            if jj == 1 && kk == 1
                
                set(gca, 'linewidth', 2, 'ylim', [1 100], 'xlim', [-100 300], ...
                    'tickdir', 'out', 'xminortick', 'on', 'yminortick', 'on')
                box off
                
                patch([-100 0 0 -100], [0 0 33.332 33.332], ...
                    'k', 'edgecolor', 'none',  'facealpha', .1)
                
                vt = plot([0 0], [0 33.322+8.333]);
                set(vt, 'color', [0 0 0 .75], 'linewidth', 2, 'linestyle', '-')
                
                patch([58 78 78 58], ...
                    [0 0 100 100], ...
                    [1 .4 0], 'edgecolor', 'none',  'facealpha', .25)
                
                vt = vline([58 78]);
                set(vt, 'color', [1 .4 0 .25], 'linewidth', 2, 'linestyle', '-')
                
                patch([mua_quarts_compiled.time_vec(1) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(end) mua_quarts_compiled.time_vec(1)], ...
                    [sig_thres_down sig_thres_down sig_thres_up sig_thres_up], ...
                    'red', 'edgecolor', 'none',  'facealpha', .25)
                
                ht = hline([sig_thres_up sig_thres_down]);
                set(ht, 'color', [1 0 0 .25], 'linewidth', 2, 'linestyle', '-')
                
            end

            plot(mua_quarts_compiled.time_vec(1:endingp), ...
                smoothdata(squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)),'movmean',5),...
                ...squeeze(mua_quarts{kk}.CS_crit_percent_nrs(jj,ii,1:endingp)), ...
                'Color', [tmrk .5], 'linewidth', 2);
            
            
            if kk == numel(mua_quarts) && jj == mua_quarts_compiled.mpi
                ht = hline(16.666);
                set(ht, 'color', [0 0 0 .5], 'linewidth', 2, 'linestyle', ':')
                set(gca, 'xlim', [-50 300])
            end
            
        end
    end
end