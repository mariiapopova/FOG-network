clear all
close all

%% Set initial conditions
%small grid
% aopt=linspace(24,27,4);
% bopt=linspace(14,18,5);
% copt=bopt;
% dopt=[0.5 0.6 0.7];
% eopt=[1 2 3 4 5];
% fopt=2.6; 

%str
% aopt=linspace(0,2.6,4);
% bopt=linspace(1,9,5);
% copt=linspace(0.044,9*0.044,5);
% dopt=linspace(1,5,5);
% eopt=1;
% fopt=1; 

%gpe
aopt=linspace(0,0.7,7);
bopt=[0.8667 1.733];
copt=[0.132 0.396];
dopt=[3 4 5];
eopt=1;
fopt=1; 

num=0;
for i=1:length(aopt)
    for j=1:length(bopt)
        for k=1:length(copt)
            for l=1:length(dopt)
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
r=randn(n,1)*2; %what is r?

%Sensorimotor cortex input to talamic cells
[Istim, timespike]=createSMCinput(tmax,dt,14,0.2); 

%BGnetwork loads Istim.mat which has all the initial conditions
save('Istim.mat','Istim','timespike','tmax','dt','v1','v2','v3','v4','v5','v6','v7','r','n');


%% Running BGnetwork.m
%For 1000msec with 10 neurons in each nucleus, each condition will take
%roughly 60sec to run.

%calculate error index
%dbstop in FOGnetwork at 294 if (ismember(i, timespikeint/dt))
aoptdum=aopt(i);
boptdum=bopt(j);
coptdum=copt(k);
doptdum=dopt(l);

h=FOGnetwork(0,0,0,aopt(i), bopt(j), copt(k), dopt(l), eopt, fopt); %healthy
pd = FOGnetwork(1,0,0,aopt(i), bopt(j), copt(k), dopt(l), eopt, fopt); %PD
%dbs=FOGnetwork(1,1,130); %PD with DBS
% h=FOGnetwork(0,0,0); %healthy
% pd = FOGnetwork(1,0,0); %PD
% dbs=FOGnetwork(1,1,130); %PD with DBS
num=num+1;
disp(num);

filename=('optimization.mat'); 
m = matfile(filename,'Writable',true);
new_row_of_values = [aoptdum, boptdum, coptdum, doptdum, eopt, fopt, h, pd];
if isprop(m, 'optimization')
    s = size(m, 'optimization');
    m.optimization(s(1)+1, :) = new_row_of_values;
else
    m.optimization = new_row_of_values;
end 

%calculate wo error index
%FOGnetwork(0,0,0); %healthy
%FOGnetwork(1,0,0); %PD
%FOGnetwork(1,1,130); %PD with DBS
            end
        end
    end
end



