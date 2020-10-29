clear all
close all

%% This file is to depict multiple comparisons results
% FR model - LC
x=categorical({'Healthy' 'PD' 'STN-DBS' 'STN+SNr-DBS'});
x=reordercats(x,{'Healthy' 'PD' 'STN-DBS' 'STN+SNr-DBS'});
lc_std=[0.9858, 0.9779, 1.3713, 0];
lc_mean=[6.2414, 7.8007, 6.7133, 0];

lc1=figure(1);
bar(x,lc_mean)
hold on
errorbar(x, lc_mean, lc_std)
hold off
ylabel('LC firing rate, Hz')
xlabel('Network states')
%saveas(lc1, "FR_lc.png")

% FR model - PPN
ppn_std=[1.3660, 1.6025, 1.8463, 1.1721];
ppn_mean=[16.5050, 21.0049, 20.4012, 11.9176];

ppn1=figure(2);
bar(x, ppn_mean)
hold on
errorbar(x, ppn_mean, ppn_std)
hold off
ylabel('PPN firing rate, Hz')
xlabel('Network states')
%saveas(ppn1, "FR_ppn.png")

%% SER model
% PPN
load('PPNactivations.mat')
ppn_pd = PPNactivations(:,1)/10000; %change with pther times
ppn_pd = squeeze(ppn_pd);
ppn_pd_mean = mean(ppn_pd);
ppn_pd_std = std(ppn_pd);

ppn_pd_stn = PPNactivations(:,2)/10000; %or 2? %change with other times
ppn_pd_stn = squeeze(ppn_pd_stn);
ppn_pd_stn_mean = mean(ppn_pd_stn);
ppn_pd_stn_std = std(ppn_pd_stn);

ppn_pd_snr = PPNactivations(:,4)/10000; %or 4? %change with other times
ppn_pd_snr = squeeze(ppn_pd_snr);
ppn_pd_snr_mean = mean(ppn_pd_snr);
ppn_pd_snr_std = std(ppn_pd_snr);

x1=categorical({'PD' 'STN-DBS' 'STN+SNr-DBS'});
x1=reordercats(x1,{'PD' 'STN-DBS' 'STN+SNr-DBS'});
ppn_mean2 = [ppn_pd_mean, ppn_pd_stn_mean, ppn_pd_snr_mean];
ppn_std2 = [ppn_pd_std, ppn_pd_stn_std, ppn_pd_snr_std];

ppn2=figure(3);
bar(x1,ppn_mean2)
hold on
errorbar(x1, ppn_mean2, ppn_std2)
hold off
ylabel('PPN average activation, %')
xlabel('Network states')
%saveas(ppn2, "SER_ppn.png")

% LC
load('LCactivations.mat')
lc_pd = LCactivations(:,1)/10000; %change with pther times
lc_pd = squeeze(lc_pd);
lc_pd_mean = mean(lc_pd);
lc_pd_std = std(lc_pd);

lc_pd_stn = LCactivations(:,2)/10000; %or 2? %change with other times
lc_pd_stn = squeeze(lc_pd_stn);
lc_pd_stn_mean = mean(lc_pd_stn);
lc_pd_stn_std = std(lc_pd_stn);

lc_pd_snr = LCactivations(:,4)/10000; %or 4? %change with other times
lc_pd_snr = squeeze(lc_pd_snr);
lc_pd_snr_mean = mean(lc_pd_snr);
lc_pd_snr_std = std(lc_pd_snr);

lc_mean2 = [lc_pd_mean, lc_pd_stn_mean, lc_pd_snr_mean];
lc_std2 = [lc_pd_std, lc_pd_stn_std, lc_pd_snr_std];

lc2=figure(4);
bar(x1,lc_mean2)
hold on
errorbar(x1, lc_mean2, lc_std2)
hold off
ylabel('LC average activation, %')
xlabel('Network states')
%saveas(lc2, "SER_ppn.png")