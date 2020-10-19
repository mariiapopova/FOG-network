clear all
close all

%% Set initial conditions

num=0;
iter=100;
for i=1:iter

%time variables
tmax=1000; %maximum time (ms)
dt=0.01; %timestep (ms)
t=0:dt:tmax; %time vector
n=10; %number of neurons in each nucleus (TH, STN, GPe, GPi)

%initial membrane voltages for all cells
v1=-62+randn(n,1)*5;
v2=-62+randn(n,1)*5;
v3=-62+randn(n,1)*5;
v4=-62+randn(n,1)*5;
v5=-62+randn(n,1)*5; %for PPN
v6=-62+randn(n,1)*5; %for SNR
v7=-63.8+randn(n,1)*5; %for striatum
v8=-62+randn(n,1)*5; %for PRF
v9=-62+randn(n,1)*5; %for CNF
v10=-62+randn(n,1)*5; %for LC
r=randn(n,1)*2; %what is r?

%Sensorimotor cortex input to talamic cells
[Istim, timespike]=createSMCinput(tmax,dt,14,0.2); 

%BGnetwork loads Istim.mat which has all the initial conditions
save('Istim.mat','Istim','timespike','tmax','dt','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','r','n');


%% Running BGnetwork.m
%For 1000msec with 10 neurons in each nucleus, each condition will take
%roughly 60sec to run.

%calculate error index
[h,fr1h,fr2h,fr3h,fr4h,fr5h,fr6h,fr7h,fr8h,fr9h,fr10h]=FOGnetwork(0,0,0); %healthy
[pd,fr1p,fr2p,fr3p,fr4p,fr5p,fr6p,fr7p,fr8p,fr9p,fr10p] = FOGnetwork(1,0,0); %PD
[dbs,fr1d,fr2d,fr3d,fr4d,fr5d,fr6d,fr7d,fr8d,fr9d,fr10d]=FOGnetwork(1,1,130); %PD with DBS

num=num+1;
disp(num);

filename=('mruns1410stn.mat'); 
m = matfile(filename,'Writable',true);
new_row_of_values = [fr1h, fr2h, fr3h, fr4h, fr5h, fr6h, fr7h,fr8h,fr9h,fr10h, h, fr1p, fr2p, fr3p, fr4p, fr5p, fr6p, fr7p,fr8p,fr9p,fr10p, pd, fr1d, fr2d, fr3d, fr4d, fr5d, fr6d, fr7d,fr8d,fr9d,fr10d, dbs];
if isprop(m, 'mruns1410stn')
    s = size(m, 'mruns1410stn');
    m.mruns1410stn(s(1)+1, :) = new_row_of_values;
else
    m.mruns1410stn = new_row_of_values;
end 

end
