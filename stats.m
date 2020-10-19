clear all
close all

%% Calculate statistics fpr a dataset
load('mruns1410snr.mat');
mruns=mruns1410snr;

%healthy data
sn_h=mruns(:,1);
ge_h=mruns(:,2);
gi_h=mruns(:,3);
ppn_h=mruns(:,4);
th_h=mruns(:,5);
str_h=mruns(:,6);
snr_h=mruns(:,7);
prf_h=mruns(:,8);
cnf_h=mruns(:,9);
lc_h=mruns(:,10);
h=mruns(:,11);

%pd data
sn_p=mruns(:,12);
ge_p=mruns(:,13);
gi_p=mruns(:,14);
ppn_p=mruns(:,15);
th_p=mruns(:,16);
str_p=mruns(:,17);
snr_p=mruns(:,18);
prf_p=mruns(:,19);
cnf_p=mruns(:,20);
lc_p=mruns(:,21);
p=mruns(:,22);

%dbs data
sn_d=mruns(:,23);
ge_d=mruns(:,24);
gi_d=mruns(:,25);
ppn_d=mruns(:,26);
th_d=mruns(:,27);
str_d=mruns(:,28);
snr_d=mruns(:,29);
prf_d=mruns(:,30);
cnf_d=mruns(:,31);
lc_d=mruns(:,32);
d=mruns(:,33);

%Calculate the means
sn_h_m=mean(sn_h);
ge_h_m=mean(ge_h);
gi_h_m=mean(gi_h);
ppn_h_m=mean(ppn_h);
th_h_m=mean(th_h);
str_h_m=mean(str_h);
snr_h_m=mean(snr_h);
prf_h_m=mean(prf_h);
cnf_h_m=mean(cnf_h);
lc_h_m=mean(lc_h);
h_m=mean(h);

sn_p_m=mean(sn_p);
ge_p_m=mean(ge_p);
gi_p_m=mean(gi_p);
ppn_p_m=mean(ppn_p);
th_p_m=mean(th_p);
str_p_m=mean(str_p);
snr_p_m=mean(snr_p);
prf_p_m=mean(prf_p);
cnf_p_m=mean(cnf_p);
lc_p_m=mean(lc_p);
p_m=mean(p);

sn_d_m=mean(sn_d);
ge_d_m=mean(ge_d);
gi_d_m=mean(gi_d);
ppn_d_m=mean(ppn_d);
th_d_m=mean(th_d);
str_d_m=mean(str_d);
snr_d_m=mean(snr_d);
prf_d_m=mean(prf_d);
cnf_d_m=mean(cnf_d);
lc_d_m=mean(lc_d);
d_m=mean(d);

%Calculate stds
sn_h_s=std(sn_h);
ge_h_s=std(ge_h);
gi_h_s=std(gi_h);
ppn_h_s=std(ppn_h);
th_h_s=std(th_h);
str_h_s=std(str_h);
snr_h_s=std(snr_h);
prf_h_s=std(prf_h);
cnf_h_s=std(cnf_h);
lc_h_s=std(lc_h);
h_s=std(h);

sn_p_s=std(sn_p);
ge_p_s=std(ge_p);
gi_p_s=std(gi_p);
ppn_p_s=std(ppn_p);
th_p_s=std(th_p);
str_p_s=std(str_p);
snr_p_s=std(snr_p);
prf_p_s=std(prf_p);
cnf_p_s=std(cnf_p);
lc_p_s=std(lc_p);
p_s=std(p);

sn_d_s=std(sn_d);
ge_d_s=std(ge_d);
gi_d_s=std(gi_d);
ppn_d_s=std(ppn_d);
th_d_s=std(th_d);
str_d_s=std(str_d);
snr_d_s=std(snr_d);
prf_d_s=std(prf_d);
cnf_d_s=std(cnf_d);
lc_d_s=std(lc_d);
d_s=std(d);

%bar plot
x=categorical({'healthy','pd','dbs'});
x=reordercats(x,{'healthy','pd','dbs'});

figure;
subplot(2,6,1);
sn=[sn_h_m, sn_p_m, sn_d_m];
bar(x,sn);
hold on
errorbar(x,sn,[sn_h_s, sn_p_s, sn_d_s],'.')
hold off
title('STN')

subplot(2,6,2);
ge=[ge_h_m, ge_p_m, ge_d_m];
bar(x,ge)
hold on
errorbar(x,ge,[ge_h_s, ge_p_s, ge_d_s],'.')
hold off
title('GPe')

subplot(2,6,3);
gi=[gi_h_m, gi_p_m, gi_d_m];
bar(x,gi)
hold on
errorbar(x,gi,[gi_h_s, gi_p_s, gi_d_s],'.')
hold off
title('GPi')

subplot(2,6,4);
str=[str_h_m, str_p_m, str_d_m];
bar(x,str)
hold on
errorbar(x,str,[str_h_s, str_p_s, str_d_s],'.')
hold off
title('Striatum')

subplot(2,6,5);
th=[th_h_m, th_p_m, th_d_m];
bar(x,th)
hold on
errorbar(x,th,[th_h_s, th_p_s, th_d_s],'.')
hold off
title('Thalamus')

subplot(2,6,6);
ppn=[ppn_h_m, ppn_p_m, ppn_d_m];
bar(x,ppn)
hold on
errorbar(x,ppn,[ppn_h_s, ppn_p_s, ppn_d_s],'.')
hold off
title('PPN')

subplot(2,6,7);
snr=[snr_h_m, snr_p_m, snr_d_m];
bar(x,snr)
hold on
errorbar(x,snr,[snr_h_s, snr_p_s, snr_d_s],'.')
hold off
title('SNr')

subplot(2,6,8);
prf=[prf_h_m, prf_p_m, prf_d_m];
bar(x,prf)
hold on
errorbar(x,prf,[prf_h_s, prf_p_s, prf_d_s],'.')
hold off
title('PRF')

subplot(2,6,9);
cnf=[cnf_h_m, cnf_p_m, cnf_d_m];
bar(x,cnf)
hold on
errorbar(x,cnf,[cnf_h_s, cnf_p_s, cnf_d_s],'.')
hold off
title('CNF')

subplot(2,6,10);
lc=[lc_h_m, lc_p_m, lc_d_m];
bar(x,lc)
hold on
errorbar(x,lc,[lc_h_s, lc_p_s, lc_d_s],'.')
hold off
title('LC')

subplot(2,6,11);
gen=[h_m, p_m, d_m];
bar(x,gen)
hold on
errorbar(x,gen,[h_s, p_s, d_s],'.')
hold off
title('Error Index')

suptitle(['Firing rates and EI in FOG network']);
disp(d_m);
disp(d_s);