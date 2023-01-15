load("Z:\_DATA\MANUSCRIPTS\WesterbergEtAl_2022_NCOMMS_Feedforward\Revision_Data\mua_targ_cor_dist_cor_compiled_for_PRAv3.mat")

targ_resp = reli_dat_pop_targ;
dist_resp = reli_dat_pop_dist;
targ_rt = rt_targ;
dist_rt = rt_dist;
targ_rnk = rtrnk_targ;
dist_rnk = rtrnk_dist;

itts = 1;

maxvins = 25;
maxpopsize = 250;
maxbin  = 151;
a_sum = nan(1000, maxvins, maxbin, itts);
a_rtm = nan(1000, maxvins, maxbin, itts);

reli_dat_pop_targ = targ_resp;
rt_targ = targ_rt;
rtrnk_targ = targ_rnk;

for iii = 1 : itts
for ii = 1:maxbin %407-25
for j = 250
    kctr = 0;
    for k = 0:(100/maxvins):99.99999999999
        kctr = kctr + 1;

        prct = [k k+(100/maxvins)];
        tdat = nanmean(reli_dat_pop_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2), 168:168+ii), 2);
        %tdat = nanmean(reli_dat_pop_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2), ii:ii+24), 2);
        trt = rt_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2));
        tsamp = datasample(1:numel(trt), j * 1000, 'Replace', true);

        a_sum(:,kctr,ii, iii) = ...
            nansum( reshape(tdat(tsamp), 1000, j), 2 );

        a_rtm(:,kctr,ii, iii) = ...
            nanmedian( reshape(trt(tsamp), 1000, j), 2 );

        clear prct tdat trt tsamp
    end
end
end
end

a_sum = a_sum(:,1:end-1,:,:);
a_rtm = a_rtm(:,1:end-1,:,:);

for iii = 1:itts
for i = 1 : size(a_sum, 3)
    try
        t1 = a_sum(:,:,i,iii);
        t2 = a_rtm(:,:,i,iii);
        [fr{i,iii}, gof{i,iii}] = power_fit(t1(:), t2(:));
    end
end
end

for iii = 1:itts
for i = 1 : size(fr,1) %maxbin
    try
        r2(i,iii) = gof{i}.adjrsquare;
        tt = coeffvalues(fr{i,iii});
        b(i,iii) = tt(2); clear tt;
    catch
        r2(i,iii) = NaN;
        b(i,iii) = NaN;
    end
end
end

% reli_dat_pop_targ = targ_resp(depth_pop_targ == 2, :);
% rt_targ = targ_rt(depth_pop_targ == 2);
% rtrnk_targ = targ_rnk(depth_pop_targ == 2);
% 
% for j = 250
%     kctr = 0;
%     for k = 0:(100/maxvins):99.99999999999
%         kctr = kctr + 1;
% 
%         prct = [k k+(100/maxvins)];
%         tdat = nanmean(reli_dat_pop_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2), 168:188), 2);
%         trt = rt_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2));
%         tsamp = datasample(1:numel(trt), j * 1000, 'Replace', true);
% 
%         g_sum(:,kctr) = ...
%             nansum( reshape(tdat(tsamp), 1000, j), 2 );
% 
%         g_rtm(:,kctr) = ...
%             nanmedian( reshape(trt(tsamp), 1000, j), 2 );
% 
%         clear prct tdat trt tsamp
%     end
% end
% 
% reli_dat_pop_targ = targ_resp(depth_pop_targ == 1, :);
% rt_targ = targ_rt(depth_pop_targ == 1);
% rtrnk_targ = targ_rnk(depth_pop_targ == 1);
% 
% for j = 250
%     kctr = 0;
%     for k = 0:(100/maxvins):99.99999999999
%         kctr = kctr + 1;
% 
%         prct = [k k+(100/maxvins)];
%         tdat = nanmean(reli_dat_pop_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2), 168:188), 2);
%         trt = rt_targ(rtrnk_targ>prct(1) & rtrnk_targ<=prct(2));
%         tsamp = datasample(1:numel(trt), j * 1000, 'Replace', true);
% 
%         s_sum(:,kctr) = ...
%             nansum( reshape(tdat(tsamp), 1000, j), 2 );
% 
%         s_rtm(:,kctr) = ...
%             nanmedian( reshape(trt(tsamp), 1000, j), 2 );
% 
%         clear prct tdat trt tsamp
%     end
% end
