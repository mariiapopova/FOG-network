function DBS=dbssyn(f,tmax,dt)
% transmission + synaptic delay: td 
td=0; %2 ms for trasmission and .5 ms for synaptic delay prev-2
ti=0; tf=tmax+td;%%in mili seconds 
t=ti:dt:tf;

% DBS input
fdbs=f;
%in mili seconds
dbsi=(dt)/dt; dbsf=tmax/dt;

% I Kernel time constant
taus=1.7;  %For excitatory synapse (Markram)
%taus=8.3;  %For inhibitory synapse (Markram)

% transmission + synaptic delay: td 
td=td/dt; %convert to simulation step scale

%input spike train
sp=zeros(1,length(t));

% Synapse parameters % Each column 1,2,3 means F,D,P respectively and each row means
% Excitatory and inhibitory synapse (1: excitatory, 2: inhibitory)
% In this study we just used the first row, excitatory synapses. (Markram)
tauf=[670,17,326; 376,21,62];
taud=[138,671,329; 45,706,144];
U=[.09,.5,.29; .016,.25,.32];
A=[.00025,.00025,.00025; .00025,.00025,.00025]; %(250pA Tsodyks)
%A=[1,1,1; 1,1,1];
n=100; A=n*A;  % change the strength of A (order of magnitude of totall number of synapses) 10 to 10
ie=ones(1,2);

fid=2.5; %synaptic fidelity 
we=fid*200; wi=0; 
% Percentage of excitatory and inhibitory synapses:
ne=we*[45,38,17]; %original: 45,38,17
% ne=zeros(1,3);
% for 1 synapse n1=1 and so forth (approximately giving 2 pA exc. current)
% ne=10;  % for 10 synapses (approximately giving 20 pA exc. current)
% ne=100; % for 100 synapses (approximately giving 200 pA exc. current)
% ne=1000;% for 1000 synapses (approximately giving 2 nA exc. current)
% ni=wi*[13,10,6];     % for 1 synapse (approximately giving 10 pA inhibitory current)
ni=wi*[8,76,16]; %ne=ni;
% ni=zeros(1,3);
% ni=10;  % for 10 synapses (approximately giving 100 pA inh. current)
% ni=100; % for 100 synapses (approximately giving 1 nA inh. current)
% ni=1000;% for 1000 synapses (approximately giving 10 nA inh. current)
A=[ne.*A(1,:);ni.*A(2,:)];

% Compute EPSC
u=zeros(1,length(t));
x=ones(1,length(t));
I=zeros(1,length(t));
PSC=zeros(length(ie),length(A),length(t));

% Compute neuron firing pattern with and without synaptic input:
DBS=zeros(1,length(t));

for q=1:length(ie)
    if q==1 
        w=1;
    else
        w=-1;
    end
for p=1:length(A)
        T=round((1000/fdbs)/dt);
        dbs=dbsi:T:dbsf;
        ts=dbs;      %uncomment for DBS only
        sp(ts)=1/dt;
            for i=td+1:length(t)-1  
                u(i+1) = u(i) + dt*(-(u(i)/tauf(q,p))+U(q,p)*(1-u(i))*sp(i-td)); 
                x(i+1) = x(i) + dt*((1/taud(q,p))*(1-x(i)) - u(i+1)*x(i)*sp(i-td));
                I(i+1) = I(i) + dt*((-1/taus)*I(i) + A(q,p)*u(i+1)*x(i)*sp(i-td));
            end
            %replace I with Iwo for no depletion
        PSC(q,p,:)= w*I; 
    end
end

PSC_exc=sum(PSC(1,:,:),2);
PSC_inh=sum(PSC(2,:,:),2);
PSC_all=PSC_exc+PSC_inh;
DBS=squeeze(PSC_all(1,1,:));


