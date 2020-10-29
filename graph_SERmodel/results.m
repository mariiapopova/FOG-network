%clear all
%close all
%% Neuronal state in PD case - SER simulation
%initials
load('FOG_directed.mat') %directed FOG network 
T = 100; %time
NN = 1000; %number of neurons
PD = [0,0,0,0,0,0,1,0,0,1]; %initial state of the regions for PD state 
                               %in order (LC,PRF,CNF,PPN,SNr,STN,GPi,GPe,Str,Ctx)                               
adj_pd = FOG_directed.'; %need to transpose as it was from python %adjacency matrix for pd patients

res_pd = SERmodel_multneuro(NN, adj_pd, T, PD, 0); %results for PD state directed graph/no dbs
%% calculate statistics for LC in PD state
lc_pd = squeeze(res_pd(:,4,:)); %take only LC part
lc_activ = (lc_pd==1); %take all neurons with 1 in all time
%lc_tot_activ = sum(sum(lc_activ)); %number of spikes in LC
lc_tot_activ = sum(lc_activ(:,end));
%% change matrix to get STN-DBS
FOG_directed_stn = FOG_directed;
FOG_directed_stn(6,:) = [0 0 0 1.5 1.5 0 1.5 1.5 0 0]; %all connections from stn 1.5
adj_stn = FOG_directed_stn.';
STN = [0,0,0,0,0,1,1,0,0,1]; %stn is active
res_stn = SERmodel_multneuro(NN, adj_stn, T, STN, 0); %results for STN state directed graph
%% calculate statistics for LC in STN-DBS state
lc_stn = squeeze(res_stn(:,4,:)); %take only LC part
lc_stn_activ = (lc_stn==1); %take all neurons with 1 in all time
%lc_stn_tot_activ = sum(sum(lc_stn_activ)); %number of spikes in LC
lc_stn_tot_activ = sum(lc_stn_activ(:,end));

%%%% try with stn always on
res_stn_on = SERmodel_multneuro(NN, adj_stn, T, STN, 1); %results for STN state directed graph
lc_stn_on = squeeze(res_stn_on(:,4,:)); %take only LC part
lc_stn_on_activ = (lc_stn_on==1); %take all neurons with 1 in all time
%lc_stn_on_tot_activ = sum(sum(lc_stn_on_activ)); %number of spikes in LC
lc_stn_on_tot_activ = sum(lc_stn_on_activ(:,end));
%% change matrix to get STN+SNr-DBS
FOG_directed_snr = FOG_directed_stn;
FOG_directed_snr(5,:) = [0 0 -1.5 -1.5 0 0 0 0 0 0]; %all connections from snr 1.5
adj_snr = FOG_directed_snr.';
SNR = [0,0,0,0,1,1,1,0,0,1]; %snr is active
res_snr = SERmodel_multneuro(NN, adj_snr, T, SNR, 0); %results for SNR state directed graph
%% calculate statistics for LC in STN-DBS state
lc_snr = squeeze(res_snr(:,4,:)); %take only LC part
lc_snr_activ = (lc_snr==1); %take all neurons with 1 in all time
%lc_snr_tot_activ = sum(sum(lc_snr_activ)); %number of spikes in LC
lc_snr_tot_activ = sum(lc_snr_activ(:,end));

%%%% try with stn always on
res_snr_on = SERmodel_multneuro(NN, adj_snr, T, SNR, 2); %results for STN state directed graph
lc_snr_on = squeeze(res_snr_on(:,4,:)); %take only LC part
lc_snr_on_activ = (lc_snr_on==1); %take all neurons with 1 in all time
%lc_snr_on_tot_activ = sum(sum(lc_snr_on_activ)); %number of spikes in LC
lc_snr_on_tot_activ=sum(lc_snr_on_activ(:,end));