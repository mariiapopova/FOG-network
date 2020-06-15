function [Sc, GN] = FOGnetwork(pd,stim,freq)

%Usage: GN = FOGnetwork(pd,stim,freq)
%
%Example: error_index=BG(1,1,130);
%
%Variables:
%pd - Variable to determine whether network is under the healthy or 
%Parkinsonian condition. For healthy, pd = 0, for Parkinson's, pd = 1.
%stim - Variable to determine whether deep brain stimulation is on.
%If DBS is off, stim = 0. If DBS is on, stim = 1.
%freq - Determines the frequency of stimulation, in Hz.
%
%
%Author: Mariia Popova, UKE; based on Rosa So, Duke University 
%Updated 9/06/2020

load('Istim.mat') %loads initial conditions %what is r?
addpath('gating')
%change the synapse to simplistic 1 to 1 structure!!!!
%Membrane parameters
Cm=1; 
%Csnr=0.1; %for SNr membrane
%In order of Th,STN,GP,PPN,SNr or Th,STN,GPe,GPi,PPN,SNr
gl=[0.05 2.25 0.1]; El=[-70 -60 -65];
gna=[3 37 120 30]; Ena=[50 55 55 45]; 
gk=[5 45 30 3.2]; Ek=[-75 -80 -80 -95];
gt=[5 0.5 0.5]; Et=0;
gca=[0 2 0.15]; Eca=[0 140 120];
gahp=[0 20 10]; %eahp and ek are the same excluding th
gnal5=0.0207; %na leakage ppn
gkl5=0.05; %k leakage ppn
gcort=0.15; Ecort=0; %cortex par for ppn
ghyp5=0.4; Ehyp5=-43; %hyperpolarization-activated current ppn
gnap5=1.9; %persistent sodium current ppn - chosen to have 8 Hz in rest
Pca=1*10^(-3); z=2; F=96490; R=8.314; T=309.5; Cao=2; Cai=2.4e-4; %alike destexhe ghk for ppn
Bcort=1;%ms^-1?

k1=[0 15 10]; %dissociation const of ahp current 
kca=[0 22.5 15]; %calcium pump rate constant
%A=[0 3 2 2 3]; %params for synapses %check with Lourens?
%B=[0 0.1 0.04 0.04 0.1]; %params for synapses %check with Lourens?
%the=[0 30 20 20 30]; %params for synapses %why so many if need only for gpe
A=[0 1 2 2 1]; %which ones to take??
B=[0 5 0.04 0.1 5]; 
the=[0 30 20 20 30];
%A=[0 3 2 2 3]; %which ones to take??
%B=[0 0.1 0.04 0.04 0.1]; 
%the=[0 30 20 20 30];

%%Synapse parameters
%In order of Igesn,Isnge,Igege,Isngi,Igegi,Igith 
gsyn = [1 0.3 1 0.3 1 .08]; Esyn = [-85 0 -85 0 -85 -85]; %alike in Rubin gsyn and in So Esyn
%gsynppn = [0.15 0.1 0.15]; Esynppn = [0 -95 0]; %alike in Lourens
gsynppn = [0.7 0.06 0.7 0.06]; Esynppn = [0 -85 0 0]; %works smh can be 0.08
%gsynppn = [0 0.06 0 0.06]; Esynppn = [0 -85 0 0]; 
tau=5; gpeak1=0.3; gpeak=0.43; %parameters for second-order alpha synapse

%time step
t=0:dt:tmax;

%%Setting initial matrices
%n - number of neurons in each population
vth=zeros(n,length(t)); %thalamic membrane voltage
vsn=zeros(n,length(t)); %STN membrane voltage
vge=zeros(n,length(t)); %GPe membrane voltage
vgi=zeros(n,length(t)); %GPi membrane voltage
vppn=zeros(n,length(t)); %PPN membrane voltage
vsnr=zeros(n,length(t));%SNr membrane voltage
Z4=zeros(n,1); %for 2 order alpha-synapse gpi-th current
S4=zeros(n,1); %for alpha-synapse gpi-th current
S3=zeros(n,1); %for alpha-synapse gesn current
S31=zeros(n,1);%for dummy gesn current
S2=zeros(n,1); %for alpha-synapse snge current
Z2=zeros(n,1); %for 2 order alpha-synapse sn current
S21=zeros(n,1); %for dummy snge current
S32=zeros(n,1); %for dummy gege current
S5=zeros(n,1); %for alpha-synapse stn-ppn current
S51=zeros(n,1); %for dummy snppn current
S52=zeros(n,1); %for dummy snppn current
S53=zeros(n,1); %for dummy snppn current
S54=zeros(n,1); %for dummy snppn current
S6=zeros(n,1); %for alpha-synapse gpi-ppn current
S61=zeros(n,1); %for dummy gippn current
S62=zeros(n,1); %for dummy gippn current
S63=zeros(n,1); %for dummy gippn current
S64=zeros(n,1); %for dummy gippn current
S7=zeros(n,1); %for alpha-synapse ppn-stn current
S71=zeros(n,1); %for dummy ppnsn current
S72=zeros(n,1); %for dummy ppnsn current
S73=zeros(n,1); %for dummy ppnsn current
S8=zeros(n,1); %for alpha-synapse ppn-gpi current
S81=zeros(n,1); %for dummy ppngi current
S82=zeros(n,1); %for dummy ppngi current
S83=zeros(n,1); %for dummy ppngi current
Sc=zeros(n,length(t)); %for cortex-ppn synapse

%%with or without DBS
Idbs=createdbs(freq,tmax,dt); %creating DBS train with frequency freq
if ~stim; Idbs=zeros(1,length(t)); end

%%initial conditions 
vth(:,1)=v1;
vsn(:,1)=v2;
vge(:,1)=v3;
vgi(:,1)=v4;
vppn(:,1)=v5; %for PPN
%vsnr(:,1)=v6; %for SNr

%helper variables for gating and synapse params - starting parameters
R2=stn_rinf(vsn(:,1)); %r for stn
H1=th_hinf(vth(:,1)); %h for th
R1=th_rinf(vth(:,1)); %r for th
N2=stn_ninf(vsn(:,1));%n for stn
H2=stn_hinf(vsn(:,1));%h for stn
C2=stn_cinf(vsn(:,1));%c for stn
CA2=0.1; %intracellular concentration of Ca2+ in muM for stn
CA3=CA2; %for gpe
CA4=CA2; %for gpi
N3=gpe_ninf(vge(:,1));%n for gpe
H3=gpe_hinf(vge(:,1));%h for gpe
R3=gpe_rinf(vge(:,1));%r for gpe
N4=gpe_ninf(vgi(:,1));%n for gpi
H4=gpe_hinf(vgi(:,1));%h for gpi
R4=gpe_rinf(vgi(:,1));%r for gpi
M5=ppn_minf(vppn(:,1));%m for ppn na
H5=ppn_hinf(vppn(:,1));%h for ppn na
Mk5=ppn_mkinf(vppn(:,1));%m for ppn k
Mh5=ppn_mhinf(vppn(:,1));%m for ppn hyp
Hnap5=ppn_hnapinf(vppn(:,1));%h for ppn nap
Mt5=ppn_mtinf(vppn(:,1));%m for ppn low-tresh
Ht5=ppn_htinf(vppn(:,1));%h for ppn low-tresh
%M6=snr_minf(vsnr(:,1));%m for snr
%H6=snr_hinf(vsnr(:,1));%h for snr
%SNA6=snr_sinf(vsnr(:,1));%s for snr for na

timespikeint = int64(timespike); %index when 1 for ppn
%%Time loop
for i=2:length(t)
    if ismember(i, timespikeint)
        Sc(:,i-1) = 1;
    end
    V1=vth(:,i-1);    V2=vsn(:,i-1);     V3=vge(:,i-1);    V4=vgi(:,i-1); V5=vppn(:,i-1); V6=vsnr(:,i-1); %previous values
    % Synapse parameters 
    S21(2:n)=S2(1:n-1);S21(1)=S2(n); %dummy synapse for snge current as there is 1 stn to 2 ge
    S31(1:n-1)=S3(2:n);S31(n)=S3(1); %dummy synapse for gesn current as there is 1 ge to 2 stn
    S32(3:n)=S3(1:n-2);S32(1:2)=S3(n-1:n); %dummy synapse for gege current as there is 1 ge to 2 ge
    S51(3:n)=S5(1:n-2);S51(1:2)=S5(n-1:n); %dummy synapse for stnppn current
    S52(2:n)=S5(1:n-1);S52(1)=S5(n); %dummy synapse for stnppn current
    S53(6:n)=S5(1:n-5);S53(1:5)=S5(n-4:n); %dummy synapse for stnppn current
    S54(5:n)=S5(1:n-4);S54(1:4)=S5(n-3:n); %dummy synapse for stnppn current
    S61(3:n)=S6(1:n-2);S61(1:2)=S6(n-1:n); %dummy synapse for gpippn current
    S62(2:n)=S6(1:n-1);S62(1)=S6(n); %dummy synapse for gpippn current
    S63(6:n)=S6(1:n-5);S63(1:5)=S6(n-4:n); %dummy synapse for gpippn current
    S64(5:n)=S6(1:n-4);S64(1:4)=S6(n-3:n); %dummy synapse for gpippn current
    S71(2:n)=S7(1:n-1);S71(1)=S7(n); %dummy synapse for ppnsn current
    S72(5:n)=S7(1:n-4);S72(1:4)=S7(n-3:n); %dummy synapse for ppnsn current
    S73(4:n)=S7(1:n-3);S73(1:3)=S7(n-2:n); %dummy synapse for ppnsn current
    S81(2:n)=S8(1:n-1);S81(1)=S8(n); %dummy synapse for ppngi current
    S82(5:n)=S8(1:n-4);S82(1:4)=S8(n-3:n); %dummy synapse for ppnsn current
    S83(4:n)=S8(1:n-3);S83(1:3)=S8(n-2:n); %dummy synapse for ppnsn current
    
    %membrane parameters - gating variables
    m1=th_minf(V1);m2=stn_minf(V2);m3=gpe_minf(V3);m4=gpe_minf(V4); %gpe and gpi are modeled similarily
    n2=stn_ninf(V2);n3=gpe_ninf(V3);n4=gpe_ninf(V4);
    h1=th_hinf(V1);h2=stn_hinf(V2);h3=gpe_hinf(V3);h4=gpe_hinf(V4);
    p1=th_pinf(V1);
    a2=stn_ainf(V2); a3=gpe_ainf(V3);a4=gpe_ainf(V4); %for low-treshold ca
    b2=stn_binf(R2);
    s3=gpe_sinf(V3);s4=gpe_sinf(V4);
    r1=th_rinf(V1);r2=stn_rinf(V2);r3=gpe_rinf(V3);r4=gpe_rinf(V4);
    c2=stn_cinf(V2);
    m5=ppn_minf(V5); h5=ppn_hinf(V5); mk5=ppn_mkinf(V5); mh5=ppn_mhinf(V5);
    mnap5=ppn_mnapinf(V5); hnap5=ppn_hnapinf(V5); mt5=ppn_mtinf(V5); ht5=ppn_htinf(V5);
    %m6=snr_minf(V6); h6=snr_hinf(V6); sna6=snr_sinf(V6);

    %membrane parameters - time constants
    tn2=stn_taun(V2);tn3=gpe_taun(V3);tn4=gpe_taun(V4);
    th1=th_tauh(V1);th2=stn_tauh(V2);th3=gpe_tauh(V3);th4=gpe_tauh(V4);
    tr1=th_taur(V1);tr2=stn_taur(V2);tr3=30;tr4=30;
    tc2=stn_tauc(V2);
    tm5=ppn_taum(V5); th5=ppn_tauh(V5); tmk5=ppn_taumk(V5); tmh5=ppn_taumh(V5);
    thnap5=ppn_tauhnap(V5); tmt5=ppn_taumt(V5); tht5=ppn_tauht(V5);
    %tm6=snr_taum(V6);th6=snr_tauh(V6);tsna6=snr_taus(V6);
    
    %thalamic cell currents
    Il1=gl(1)*(V1-El(1));
    Ina1=gna(1)*(m1.^3).*H1.*(V1-Ena(1));
    Ik1=gk(1)*((0.75*(1-H1)).^4).*(V1-Ek(1)); %misspelled in So paper
    It1=gt(1)*(p1.^2).*R1.*(V1-Et);
    Igith=1.4*gsyn(6)*(V1-Esyn(6)).*S4; %for alpha-synapse second order kinetics
    
    %STN cell currents
    Il2=gl(2)*(V2-El(2));
    Ik2=gk(2)*(N2.^4).*(V2-Ek(2));
    Ina2=gna(2)*(m2.^3).*H2.*(V2-Ena(2));
    It2=gt(2)*(a2.^3).*(b2.^2).*(V2-Eca(2)); %misspelled in So paper
    Ica2=gca(2)*(C2.^2).*(V2-Eca(2));
    Iahp2=gahp(2)*(V2-Ek(2)).*(CA2./(CA2+k1(2))); %cause ek and eahp are the same
    Igesn=0.5*(gsyn(1)*(V2-Esyn(1)).*(S3+S31)); %first-order kinetics 1ge to 2sn
    Iappstn=33-pd*10; %to switch between parkinsonian and healthy state
    Ippnsn=0.25*(gsynppn(3)*(V2-Esynppn(3)).*(S71+S72+S7+S73)); %first-order kinetics 4ppn to stn
    
    %GPe cell currents
    Il3=gl(3)*(V3-El(3));
    Ik3=gk(3)*(N3.^4).*(V3-Ek(3));
    Ina3=gna(3)*(m3.^3).*H3.*(V3-Ena(3));
    It3=gt(3)*(a3.^3).*R3.*(V3-Eca(3)); %Eca as in Rubin and Terman
    Ica3=gca(3)*(s3.^2).*(V3-Eca(3)); %misspelled in So paper
    Iahp3=gahp(3)*(V3-Ek(3)).*(CA3./(CA3+k1(3))); %as Ek is the same with Eahp
    Isnge=0.5*(gsyn(2)*(V3-Esyn(2)).*(S2+S21)); %second-order kinetics 1sn to 2ge
    Igege=0.5*(gsyn(3)*(V3-Esyn(3)).*(S31+S32)); %first-order kinetics 1ge to 2ge
    Iappgpe=21-13*pd+r; %why need r??

    %GPi cell currents
    Il4=gl(3)*(V4-El(3));
    Ik4=gk(3)*(N4.^4).*(V4-Ek(3));
    Ina4=gna(3)*(m4.^3).*H4.*(V4-Ena(3)); %Eca as in Rubin and Terman
    It4=gt(3)*(a4.^3).*R4.*(V4-Eca(3)); %misspelled in So paper
    Ica4=gca(3)*(s4.^2).*(V4-Eca(3)); %as Ek is the same with Eahp
    Iahp4=gahp(3)*(V4-Ek(3)).*(CA4./(CA4+k1(3)));
    Isngi=0.5*(gsyn(4)*(V4-Esyn(4)).*(S2+S21)); %second-order kinetics 1sn to 2gi
    Igegi=0.5*(gsyn(5)*(V4-Esyn(5)).*(S31+S32)); %first-order kinetics 1ge to 2gi
    Iappgpi=22-pd*6;
    Ippngi=0.25*(gsynppn(4)*(V4-Esynppn(4)).*(S81+S82+S8+S83)); %first-order kinetics 4ppn to gpi
    
    %PPN cell currents
    Inal5=gnal5*(V5-Ena(4));
    Ikl5=gkl5*(V5-Ek(4));
    Ina5=gna(4)*(M5.^3).*H5.*(V5-Ena(4));
    Ik5=gk(4)*(Mk5.^4).*(V5-Ek(4));
    Ihyp5=ghyp5*(Mh5.^3).*(V5-Ehyp5);
    %alike Rubin
    Inap5=gnap5*mnap5.*Hnap5.*(V5-Ena(4)); 
    %It alike Destexhe
    zet=Pca*F*z*V5./(R*T);
    if abs(zet)>1e-4
        Gt5=(Pca*z*F)*(Cai*(-zet./(exp(-zet)-1))-Cao*(zet./(exp(zet)-1)));
    else
        Gt5=(Pca*z*F)*(Cai*(1+zet./2)-Cao*(1-zet./2));
    end
    It5=.2e-3*Mt5.^2.*Ht5.*Gt5; %alike destexhe
    Iappppn=0; %chosen from fig.5 in Lourens
    Isnppn=(gsynppn(1)*(V5-Esynppn(1)).*(S51+S52+S53+S54)); %first-order kinetics 4stn to ppn
    Igippn=(gsynppn(2)*(V5-Esynppn(2)).*(S61+S62+S63+S64)); %first-order kinetics 4gpi to ppn
    Icort=gcort*Sc(:,i-1).*(V5-Ecort);
    
    %SNr cell currents
    %Il6=gl(5)*(V6-El(5));
    %Ina6=gna(5)*(M6.^3).*H6.*SNA6*(V6-Ena(5));
    %Isk6= gsk(1)*MSK*(V6-Ek(6));
    
    %Differential Equations for cells using forward Euler method
    %thalamic
    vth(:,i)= V1+dt*(1/Cm*(-Il1-Ik1-Ina1-It1-Igith+Istim(i)));
    H1=H1+dt*((h1-H1)./th1);
    R1=R1+dt*((r1-R1)./tr1);
    Sc(:,i)=Sc(:,i-1)+dt*(-Bcort.*Sc(:,i-1));
    
    %STN
    vsn(:,i)=V2+dt*(1/Cm*(-Il2-Ik2-Ina2-It2-Ica2-Iahp2-Igesn+Iappstn+Idbs(i)-Ippnsn)); %currently STN-DBS
    N2=N2+dt*(0.75*(n2-N2)./tn2); 
    H2=H2+dt*(0.75*(h2-H2)./th2);
    R2=R2+dt*(0.2*(r2-R2)./tr2);
    CA2=CA2+dt*(3.75*10^-5*(-Ica2-It2-kca(2)*CA2));
    C2=C2+dt*(0.08*(c2-C2)./tc2); 
    %for second order second-order alpha-synapse
    a=find(vsn(:,i-1)<-10 & vsn(:,i)>-10);
    u=zeros(n,1); u(a)=gpeak/(tau*exp(-1))/dt; 
    S2=S2+dt*Z2; 
    zdot=u-2/tau*Z2-1/(tau^2)*S2;
    Z2=Z2+dt*zdot;
    %for stn-ppn synapse
    S5=S5+dt*(A(2)*(1-S5).*Hinf(V2-the(2))-B(2)*S5);
    
    %GPe
    vge(:,i)=V3+dt*(1/Cm*(-Il3-Ik3-Ina3-It3-Ica3-Iahp3-Isnge-Igege+Iappgpe));
    N3=N3+dt*(0.1*(n3-N3)./tn3); %misspelled in So paper
    H3=H3+dt*(0.05*(h3-H3)./th3); %misspelled in So paper
    R3=R3+dt*(1*(r3-R3)./tr3); %misspelled in So paper
    CA3=CA3+dt*(1*10^-4*(-Ica3-It3-kca(3)*CA3));
    S3=S3+dt*(A(3)*(1-S3).*Hinf(V3-the(3))-B(3)*S3);
    
    %GPi
    vgi(:,i)=V4+dt*(1/Cm*(-Il4-Ik4-Ina4-It4-Ica4-Iahp4-Isngi-Igegi+Iappgpi-Ippngi));
    N4=N4+dt*(0.1*(n4-N4)./tn4); %misspelled in So paper
    H4=H4+dt*(0.05*(h4-H4)./th4); %misspelled in So paper
    R4=R4+dt*(1*(r4-R4)./tr4); %misspelled in So paper
    CA4=CA4+dt*(1*10^-4*(-Ica4-It4-kca(3)*CA4));
    %for second order second-order alpha-synapse
    a=find(vgi(:,i-1)<-10 & vgi(:,i)>-10);
    u=zeros(n,1); u(a)=gpeak1/(tau*exp(-1))/dt; 
    S4=S4+dt*Z4; 
    zdot=u-2/tau*Z4-1/(tau^2)*S4;
    Z4=Z4+dt*zdot;
    %for gpi-ppn synapse
    S6=S6+dt*(A(4)*(1-S6).*Hinf(V4-the(4))-B(4)*S6);
    
    %PPN
    %vppn(:,i)=V5+dt*(1/Cm*(-Inal5-Ikl5-Ina5-Ik5-It5-Ihyp5-Inap5-Isnppn-Igippn+Iappppn+Istim(i))); %stimulation SMC current check how to deal with th and that it's like in So
    vppn(:,i)=V5+dt*(1/Cm*(-Inal5-Ikl5-Ina5-Ik5-It5-Ihyp5-Inap5+Iappppn-Isnppn-Igippn+Icort));
    H5=H5+dt*((h5-H5)./th5);
    M5=M5+dt*((m5-M5)./tm5);
    Mk5=Mk5+dt*((mk5-Mk5)./tmk5);
    Mh5=Mh5+dt*((mh5-Mh5)./tmh5);
    Mt5=Mt5+dt*((mt5-Mt5)./tmt5);
    Ht5=Ht5+dt*((ht5-Ht5)./tht5);
    Hnap5=Hnap5+dt*((hnap5-Hnap5)./thnap5); %alike rubin
    %for ppn-stn synapse
    S7=S7+dt*(A(5)*(1-S7).*Hinf(V5-the(5))-B(5)*S7);
    %current into GPi!!! for us - change GPi
    %for ppn-gpi synapse
    S8=S8+dt*(A(5)*(1-S8).*Hinf(V5-the(5))-B(5)*S8);
    
    %SNr
    %vsnr(:,i)=V6+dt*(1/Csnr*(-Il6-Ina6-Inap6-Ik6-Ica6-Isk6-Isnsnr+Iappsnr));
    %M6=M6+dt*((m6-M6)./tm6); 
    %H6=H6+dt*((h6-H6)./th6); 
    %SNA6=SNA6+dt*((sna6-SNA6)./tsna6); 
    
end

%%Calculation of error index
GN=calculateEI(t,vth,timespike,tmax); %for thalamus
%GN=calculateEI(t,vppn,timespike,tmax); %for PPN

%disp(vppn(1,:));
%figure
%plot(vppn(1,:));
%%Plots membrane potential for one cell in each nucleus
plotpotentials; 

fr=findfreq(vppn(1,:));
disp(fr);
return
