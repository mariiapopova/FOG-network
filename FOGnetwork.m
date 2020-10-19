%function GN = FOGnetwork(pd,stim,freq, aopt, bopt, copt)
%check maybe to also change gegi for pd sstate
function [GN, fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8, fr9, fr10] = FOGnetwork(pd,stim,freq)

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
%Updated 25/09/2020

load('Istim.mat') %loads initial conditions %what is r?
addpath('gating')

%Membrane parameters
Cm=1; 
%In order of Th,STN,GP,PPN,Str,PRF,LC or Th,STN,GPe,GPi,PPN,SNr,Str,PRF,LC
gl=[0.05 2.25 0.1 0.1 0.4 0.3]; El=[-70 -60 -65 -67 -70 -17];
gna=[3 37 120 30 100 120]; Ena=[50 55 55 45 50 45]; 
gk=[5 45 30 3.2 80 10 20]; Ek=[-75 -80 -80 -95 -100 -95 -72];
gt=[5 0.5 0.5]; Et=0;
gca=[0 2 0.15]; Eca=[0 140 120];
gahp=[0 20 10]; %eahp and ek are the same excluding th
gnal5=0.0207; %na leakage ppn
gkl5=0.05; %k leakage ppn
gcort=0.15; Ecort=0; %cortex par for ppn
ghyp5=0.4; Ehyp5=-43; %hyperpolarization-activated current ppn
gnap5=1.9; %persistent sodium current ppn - chosen to have 8 Hz in rest
Pca=1*10^(-3); z=2; F=96490; R=8.314; T=309.5; Cao=2; Cai=2.4e-4; %alike destexhe ghk for ppn
Bcort=1;%ms^-1
gm=1; Em=-100;%for striatum muscarinic current
%lc params
ga=47.7; Blc=0.21.*(ga/gk(7));

k1=[0 15 10]; %dissociation const of ahp current 
kca=[0 22.5 15]; %calcium pump rate constant
%synapse params alike in rubin SNr same to Gpi, PPN same to STN, Str same
%to Gpi, CNF same to PPN, PRF same to PPN - why?, LC same to PPN - why?
A=[0 5 2 2 5 2 2 5 5 5]; 
B=[0 1 0.04 0.08 1 0.08 0.08 1 1 1]; 
the=[0 30 20 20 30 20 20 30 30 30];

%%Synapse parameters
%In order of Igesn,Isnge,Igege,Isngi,Igegi,Igith 
gsyn = [1 0.3 1 0.3 1 .08]; Esyn = [-85 0 -85 0 -85 -85]; %alike in Rubin gsyn and in So Esyn
%gsynppn = [0.22 0.2 0.22 0.2 0.35]; Esynppn = [0 -85 0 0 -85]; %in order snppn gippn ppnsn ppngi snrppn can change
gsynppn = [0.2 0.18 0.2 0.18 0.18 0.18]; Esynppn = [0 -85 0 0 -85 0]; %in order snppn gippn ppnsn ppngi snrppn ppnprf
tau=5; gpeak1=0.3; gpeak=0.43; %parameters for second-order alpha synapse
gsynsnr=0.3; Esynsnr=0; %for snr synapses in order stn
gsynstr=[0.5 0.5 0.5 0.5]; ggaba=0.1; gcorstr=0.07; Esynstr=[-85 0 -85 -85 -85 -85]; tau_i=13;%parameters for striatum synapses in order gaba-rec crtx strge gestr strgi strsnr
gsynlc=[0.2 0.2 0.2]; Esynlc=[0 0 0]; %why these values?!!! - alike PPN?

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
vstr=zeros(n,length(t));%striatum membrane voltage
vprf=zeros(n,length(t));%prf membrane voltage
vcnf=zeros(n,length(t));%cnf membrane voltage
vlc=zeros(n,length(t));%LC membrane voltage

Z4=zeros(n,1); %for 2 order alpha-synapse gpi-th current
S4=zeros(n,1); %for alpha-synapse gpi-th current
S3=zeros(n,1); %for alpha-synapse gesn current
S31=zeros(n,1);%for dummy gesn current
S2=zeros(n,1); %for alpha-synapse snge current
Z2=zeros(n,1); %for 2 order alpha-synapse sn current
S21=zeros(n,1); %for dummy snge current
S32=zeros(n,1); %for dummy gege current
S5=zeros(n,1); %for alpha-synapse stn-ppn current
S6=zeros(n,1); %for alpha-synapse gpi-ppn current
S7=zeros(n,1); %for alpha-synapse ppn-stn current
S8=zeros(n,1); %for alpha-synapse ppn-gpi current
Sc=zeros(n,length(t)); %for cortex-ppn synapse
S9=zeros(n,1); %for alpha-synapse stn-snr current
S1c=zeros(n,1); %for striatum gaba-rec
%for synapses striatum-gaba-rec
all=randsample(n,n);
bll=randsample(n,n);
cll=randsample(n,n);
dll=randsample(n,n);
S10=zeros(n,1); %for alpha-synapse str-ge/gi/snr current
S11=zeros(n,1); %for alpha-synapse ge-str current
S14=zeros(n,1); %for alpha-synapse snr-ppn current
S15=zeros(n,1); %for alpha-synapse ppn-prf current
S16=zeros(n,1); %for alpha-synapse gpi-cnf current
S17=zeros(n,1); %for alpha-synapse snr-cnf current
S18=zeros(n,1); %for alpha-synapse cnf-lc current
S19=zeros(n,1); %for alpha-synapse prf-lc current
S20=zeros(n,1); %for alpha-synapse lc-prf current

%%with or without DBS
%Idbs=createdbs(freq,tmax,dt); %creating DBS train with frequency freq
Idbs=dbssyn(freq,tmax,dt); %creating DBS train with frequency freq
if ~stim; Idbs=zeros(1,length(t)); end

%%initial conditions 
vth(:,1)=v1;
vsn(:,1)=v2;
vge(:,1)=v3;
vgi(:,1)=v4;
vppn(:,1)=v5; %for PPN
vsnr(:,1)=v6; %for SNr
vstr(:,1)=v7; %for striatum D2
vprf(:,1)=v8; %for PRF
vcnf(:,1)=v9; %for CNF
vlc(:,1)=v10; %for LC

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
CA6=CA2; %for snr
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
N6=gpe_ninf(vsnr(:,1));%n for snr
H6=gpe_hinf(vsnr(:,1));%h for snr
R6=gpe_rinf(vsnr(:,1));%r for snr
N8=gpe_ninf(vprf(:,1));%n for prf
H8=gpe_hinf(vprf(:,1));%h for prf
R8=gpe_rinf(vprf(:,1));%r for prf
M9=ppn_minf(vcnf(:,1));%m for cnf na
H9=ppn_hinf(vcnf(:,1));%h for cnf na
Mk9=ppn_mkinf(vcnf(:,1));%m for cnf k
Mh9=ppn_mhinf(vcnf(:,1));%m for cnf hyp
Hnap9=ppn_hnapinf(vcnf(:,1));%h for cnf nap
Mt9=ppn_mtinf(vcnf(:,1));%m for cnf low-tresh
Ht9=ppn_htinf(vcnf(:,1));%h for cnf low-tresh
q=lc_qinf(vlc(:,1));%q for lc

%striatum gating
m7=str_alpham(vstr(:,1))./(str_alpham(vstr(:,1))+str_betam(vstr(:,1)));
h7=str_alphah(vstr(:,1))./(str_alphah(vstr(:,1))+str_betah(vstr(:,1)));
n7=str_alphan(vstr(:,1))./(str_alphan(vstr(:,1))+str_betan(vstr(:,1)));
p7=str_alphap(vstr(:,1))./(str_alphap(vstr(:,1))+str_betap(vstr(:,1)));

timespikeint = int64(timespike); %index when 1 for ppn

%%Time loop
for i=2:length(t)
    
    %condition for cortex current for ppn and striatum
    if ismember(i, timespikeint/dt)
        Sc(:,i-1) = 1;
    end
    
    V1=vth(:,i-1);    V2=vsn(:,i-1);     V3=vge(:,i-1);    V4=vgi(:,i-1);   V5=vppn(:,i-1); V6=vsnr(:,i-1); %previous values
    V7=vstr(:,i-1);   V8=vprf(:,i-1);    V9=vcnf(:,i-1);   V10=vlc(:,i-1);
    % Synapse parameters 
    S21(2:n)=S2(1:n-1);S21(1)=S2(n); %dummy synapse for snge current as there is 1 stn to 2 ge
    S31(1:n-1)=S3(2:n);S31(n)=S3(1); %dummy synapse for gesn current as there is 1 ge to 2 stn
    S32(3:n)=S3(1:n-2);S32(1:2)=S3(n-1:n); %dummy synapse for gege current as there is 1 ge to 2 ge
    S11cr=S1c(all); S12cr=S1c(bll); S13cr=S1c(cll); S14cr=S1c(dll); %dummy striatum crtx current
    
    %membrane parameters - gating variables
    m1=th_minf(V1);m2=stn_minf(V2);m3=gpe_minf(V3);m4=gpe_minf(V4);m6=gpe_minf(V6);m8=gpe_minf(V8);%gpe and gpi are modeled similarily
    n2=stn_ninf(V2);n3=gpe_ninf(V3);n4=gpe_ninf(V4);n6=gpe_ninf(V6);n8=gpe_ninf(V8);
    h1=th_hinf(V1);h2=stn_hinf(V2);h3=gpe_hinf(V3);h4=gpe_hinf(V4);h6=gpe_hinf(V6);h8=gpe_hinf(V8);
    p1=th_pinf(V1);
    a2=stn_ainf(V2); a3=gpe_ainf(V3);a4=gpe_ainf(V4);a6=gpe_ainf(V6);a8=gpe_ainf(V8); %for low-treshold ca
    b2=stn_binf(R2);
    s3=gpe_sinf(V3);s4=gpe_sinf(V4); s6=gpe_sinf(V6);
    r1=th_rinf(V1);r2=stn_rinf(V2);r3=gpe_rinf(V3);r4=gpe_rinf(V4);r6=gpe_rinf(V6);r8=gpe_rinf(V8);
    c2=stn_cinf(V2);
    m5=ppn_minf(V5); h5=ppn_hinf(V5); mk5=ppn_mkinf(V5); mh5=ppn_mhinf(V5);
    mnap5=ppn_mnapinf(V5); hnap5=ppn_hnapinf(V5); mt5=ppn_mtinf(V5); ht5=ppn_htinf(V5);
    m9=ppn_minf(V9); h9=ppn_hinf(V9); mk9=ppn_mkinf(V9); mh9=ppn_mhinf(V9);
    mnap9=ppn_mnapinf(V9); hnap9=ppn_hnapinf(V9); mt9=ppn_mtinf(V9); ht9=ppn_htinf(V9);
    qinf=lc_qinf(V10); m10=lc_minf(V10); b10=lc_binf(V10);

    %membrane parameters - time constants
    tn2=stn_taun(V2);tn3=gpe_taun(V3);tn4=gpe_taun(V4);tn6=gpe_taun(V6);tn8=gpe_taun(V8);
    th1=th_tauh(V1);th2=stn_tauh(V2);th3=gpe_tauh(V3);th4=gpe_tauh(V4);th6=gpe_tauh(V6);th8=gpe_tauh(V8);
    tr1=th_taur(V1);tr2=stn_taur(V2);tr3=30;tr4=30;tr6=30;tr8=30;
    tc2=stn_tauc(V2);
    tm5=ppn_taum(V5); th5=ppn_tauh(V5); tmk5=ppn_taumk(V5); tmh5=ppn_taumh(V5);
    thnap5=ppn_tauhnap(V5); tmt5=ppn_taumt(V5); tht5=ppn_tauht(V5);
    tm9=ppn_taum(V9); th9=ppn_tauh(V9); tmk9=ppn_taumk(V9); tmh9=ppn_taumh(V9);
    thnap9=ppn_tauhnap(V9); tmt9=ppn_taumt(V9); tht9=ppn_tauht(V9);
    tq=lc_tauq(V10); 
    
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
    Iappstn=36; 
    Ippnsn=gsynppn(3)*(V2-Esynppn(3)).*S7; %first-order kinetics ppn to stn
    
    %GPe cell currents
    Il3=gl(3)*(V3-El(3));
    Ik3=gk(3)*(N3.^4).*(V3-Ek(3));
    Ina3=gna(3)*(m3.^3).*H3.*(V3-Ena(3));
    It3=gt(3)*(a3.^3).*R3.*(V3-Eca(3)); %Eca as in Rubin and Terman
    Ica3=gca(3)*(s3.^2).*(V3-Eca(3)); %misspelled in So paper
    Iahp3=gahp(3)*(V3-Ek(3)).*(CA3./(CA3+k1(3))); %as Ek is the same with Eahp
    Isnge=0.5*(gsyn(2)*(V3-Esyn(2)).*(S2+S21)); %second-order kinetics 1sn to 2ge
    %Igege=0.5*((gsyn(3)+8*pd)*(V3-Esyn(3)).*(S31+S32));
    Igege=0.5*((gsyn(3)+1*pd)*(V3-Esyn(3)).*(S31+S32));
    Iappgpe=18;
    %str-gpe synapse
    Istrge=gsynstr(1)*(V3-Esynstr(3)).*S10; %1str to 1ge

    %GPi cell currents
    Il4=gl(3)*(V4-El(3));
    Ik4=gk(3)*(N4.^4).*(V4-Ek(3));
    Ina4=gna(3)*(m4.^3).*H4.*(V4-Ena(3)); %Eca as in Rubin and Terman
    It4=gt(3)*(a4.^3).*R4.*(V4-Eca(3)); %misspelled in So paper
    Ica4=gca(3)*(s4.^2).*(V4-Eca(3)); 
    Iahp4=gahp(3)*(V4-Ek(3)).*(CA4./(CA4+k1(3))); %as Ek is the same with Eahp
    Isngi=0.5*(gsyn(4)*(V4-Esyn(4)).*(S2+S21)); %second-order kinetics 1sn to 2gi
    Igegi=0.5*(gsyn(5)*(V4-Esyn(5)).*(S31+S32)); %first-order kinetics 1ge to 2gi
    Iappgpi=18;
    Ippngi=gsynppn(4)*(V4-Esynppn(4)).*S8; %first-order kinetics ppn to gpi
    %str-gpi synapse
    Istrgi=gsynstr(3)*(V4-Esynstr(5)).*S10; %1str to 1gi
    
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
    Isnppn=gsynppn(1)*(V5-Esynppn(1)).*S5; %first-order kinetics stn to ppn
    Igippn=gsynppn(2)*(V5-Esynppn(2)).*S6; %first-order kinetics gpi to ppn
    Icort=gcort*Sc(:,i-1).*(V5-Ecort);
    Isnrppn=gsynppn(5)*(V5-Esynppn(5)).*S14; %first-order kinetics snr to ppn
    
    %SNr cell currents - modelled as GPi
    Il6=gl(3)*(V6-El(3));
    Ik6=gk(3)*(N6.^4).*(V6-Ek(3));
    Ina6=gna(3)*(m6.^3).*H6.*(V6-Ena(3)); %Eca as in Rubin and Terman
    It6=gt(3)*(a6.^3).*R6.*(V6-Eca(3)); %misspelled in So paper
    Ica6=gca(3)*(s6.^2).*(V6-Eca(3)); 
    Iahp6=gahp(3)*(V6-Ek(3)).*(CA6./(CA6+k1(3))); %as Ek is the same with Eahp
    %Isnsnr=gsynsnr*(V6-Esynsnr).*S9; %first-order kinetics 1sn to 1snr
    Isnsnr=0.5*(gsynsnr*(V6-Esynsnr).*(S2+S21));
    %str-snr synapse
    Istrsnr=gsynstr(4)*(V6-Esynstr(6)).*S10; %1str to 1ge
    Iappsnr=0;
    
    %Striatum D2 cell currents
    Ina7=gna(5)*(m7.^3).*h7.*(V7-Ena(5));
    Ik7=gk(5)*(n7.^4).*(V7-Ek(5));
    Il7=gl(4)*(V7-El(4));
    %Im7=(2.6-2.5*pd)*gm*p7.*(V7-Em);
    Im7=(2.6-2.5*pd)*gm*p7.*(V7-Em);
    %Im7=(2.6-bopt*pd)*gm*p7.*(V7-Em);
    Igaba7=(ggaba/4)*(V7-Esynstr(1)).*(S11cr+S12cr+S13cr+S14cr); %maybe change to 3.5 for direct and indirect %recieves input from 40% remaining
    %Icorstr=(6*gcorstr-0.37*pd)*(V7-Esynstr(2)).*Sc(:,i-1); %optimized
    Icorstr=(6*gcorstr-0.3*pd)*(V7-Esynstr(2)).*Sc(:,i-1); %optimized 
    %ge-str synapse
    Igestr=gsynstr(2)*(V7-Esynstr(4)).*S11; %1ge to 1str
    Iappstr=3;
    
    %PRF cell currents - alike gpe
    Il8=gl(5)*(V8-El(5));
    Ik8=gk(6)*(N8.^4).*(V8-Ek(6));
    Ina8=gna(6)*(m8.^3).*H8.*(V8-Ena(6)); %as in Miura constants
    %Il8=gl(3)*(V8-El(3));
    %Ik8=gk(3)*(N8.^4).*(V8-Ek(3));
    %Ina8=gna(3)*(m8.^3).*H8.*(V8-Ena(3)); %as in So constants
    It8=gt(3)*(a8.^3).*R8.*(V8-Eca(3)); %Eca as in Rubin and Terman
    Icorprf=gcort*Sc(:,i-1).*(V8-Ecort); %alike ppn
    Iappprf=0;
    Ippnprf=gsynppn(6)*(V8-Esynppn(6)).*S15; %first-order kinetics ppn to gpi alike its not cholinergic
    Ilcprf=gsynlc(3)*(V8-Esynlc(3)).*S20;
    
    %CNF cell currents - alike PPN
    Inal9=gnal5*(V9-Ena(4));
    Ikl9=gkl5*(V9-Ek(4));
    Ina9=gna(4)*(M9.^3).*H9.*(V9-Ena(4));
    Ik9=gk(4)*(Mk9.^4).*(V9-Ek(4));
    Ihyp9=ghyp5*(Mh9.^3).*(V9-Ehyp5);
    %alike Rubin
    Inap9=gnap5*mnap9.*Hnap9.*(V9-Ena(4)); 
    %It alike Destexhe
    zet1=Pca*F*z*V9./(R*T);
    if abs(zet1)>1e-4
        Gt9=(Pca*z*F)*(Cai*(-zet1./(exp(-zet1)-1))-Cao*(zet1./(exp(zet1)-1)));
    else
        Gt9=(Pca*z*F)*(Cai*(1+zet1./2)-Cao*(1-zet1./2));
    end
    It9=.2e-3*Mt9.^2.*Ht9.*Gt9; %alike destexhe
    Iappcnf=0; 
    Igicnf=gsynppn(2)*(V9-Esynppn(2)).*S16; %first-order kinetics gpi to cnf alike with ppn
    Isnrcnf=gsynppn(5)*(V9-Esynppn(5)).*S17; %first-order kinetics snr to cnf
    
    %LC cell currents
    Ina10=gna(3)*(m10.^3).*(-3.*(q-Blc.*b10)+0.85).*(V10-Ena(3));
    Ik10=gk(7)*q.*(V10-Ek(7));
    Il10=gl(6).*(V10-El(6));
    Iapplc=4.5;%5
    Icnflc=gsynlc(1)*(V10-Esynlc(1)).*S18;
    Iprflc=gsynlc(2)*(V10-Esynlc(2)).*S19;
    
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
    %for stn-snr synapse
    S9=S9+dt*(A(2)*(1-S9).*Hinf(V2-the(2))-B(2)*S9);
    
    %GPe
    vge(:,i)=V3+dt*(1/Cm*(-Il3-Ik3-Ina3-It3-Ica3-Iahp3-Isnge-Igege+Iappgpe-Istrge));
    N3=N3+dt*(0.1*(n3-N3)./tn3); %misspelled in So paper
    H3=H3+dt*(0.05*(h3-H3)./th3); %misspelled in So paper
    R3=R3+dt*(1*(r3-R3)./tr3); %misspelled in So paper
    CA3=CA3+dt*(1*10^-4*(-Ica3-It3-kca(3)*CA3));
    S3=S3+dt*(A(3)*(1-S3).*Hinf(V3-the(3))-B(3)*S3);
    %ge-str synapse
    S11=S11+dt*(A(3)*(1-S11).*Hinf(V3-the(3))-B(3)*S11); 
    
    %GPi
    vgi(:,i)=V4+dt*(1/Cm*(-Il4-Ik4-Ina4-It4-Ica4-Iahp4+Iappgpi-Isngi-Igegi-Ippngi-Istrgi));
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
    %for gpi-cnf synapse
    S16=S16+dt*(A(4)*(1-S16).*Hinf(V4-the(4))-B(4)*S16);
    
    %PPN
    vppn(:,i)=V5+dt*(1/Cm*(-Inal5-Ikl5-Ina5-Ik5-It5-Ihyp5-Inap5+Iappppn-Icort-Isnppn-Igippn-Isnrppn));
    H5=H5+dt*((h5-H5)./th5);
    M5=M5+dt*((m5-M5)./tm5);
    Mk5=Mk5+dt*((mk5-Mk5)./tmk5);
    Mh5=Mh5+dt*((mh5-Mh5)./tmh5);
    Mt5=Mt5+dt*((mt5-Mt5)./tmt5);
    Ht5=Ht5+dt*((ht5-Ht5)./tht5);
    Hnap5=Hnap5+dt*((hnap5-Hnap5)./thnap5); %alike rubin
    %for ppn-stn synapse
    S7=S7+dt*(A(5)*(1-S7).*Hinf(V5-the(5))-B(5)*S7);
    %for ppn-gpi synapse
    S8=S8+dt*(A(5)*(1-S8).*Hinf(V5-the(5))-B(5)*S8);
    %for ppn-prf synapse
    S15=S15+dt*(A(5)*(1-S15).*Hinf(V5-the(5))-B(5)*S15);
    
    %SNr
    %vsnr(:,i)=V6+dt*(1/Cm*(-Il6-Ik6-Ina6-It6-Ica6-Iahp6-Isnsnr-Istrsnr+Idbs(i)+Iappsnr)); 
    vsnr(:,i)=V6+dt*(1/Cm*(-Il6-Ik6-Ina6-It6-Ica6-Iahp6+Iappsnr-Isnsnr-Istrsnr));
    N6=N6+dt*(0.1*(n6-N6)./tn6); %misspelled in So paper
    H6=H6+dt*(0.05*(h6-H6)./th6); %misspelled in So paper
    R6=R6+dt*(1*(r6-R6)./tr6); %misspelled in So paper
    CA6=CA6+dt*(1*10^-4*(-Ica6-It6-kca(3)*CA6));
    S14=S14+dt*(A(6)*(1-S14).*Hinf(V6-the(6))-B(6)*S14); %snr-ppn synapse
    S17=S17+dt*(A(6)*(1-S17).*Hinf(V6-the(6))-B(6)*S17); %snr-ppn synapse
    
    %Striatum D2
    vstr(:,i)=V7+(dt/Cm)*(-Ina7-Ik7-Il7-Im7-Igaba7-Icorstr+Iappstr-Igestr);
    m7=m7+dt*(str_alpham(V7).*(1-m7)-str_betam(V7).*m7);
    h7=h7+dt*(str_alphah(V7).*(1-h7)-str_betah(V7).*h7);
    n7=n7+dt*(str_alphan(V7).*(1-n7)-str_betan(V7).*n7);
    p7=p7+dt*(str_alphap(V7).*(1-p7)-str_betap(V7).*p7);
    S1c=S1c+dt*((str_Ggaba(V7).*(1-S1c))-(S1c/tau_i));
    %for str-gpe/gpi/str synapse
    S10=S10+dt*(A(7)*(1-S10).*Hinf(V7-the(7))-B(7)*S10); 
    
    %PRF - HH with T-current taken from GPe
    vprf(:,i)=V8+dt*(1/Cm*(-Il8-Ik8-Ina8-It8+Iappprf-Icorprf-Ippnprf-Ilcprf));
    N8=N8+dt*(0.1*(n8-N8)./tn8); 
    H8=H8+dt*(0.05*(h8-H8)./th8); 
    R8=R8+dt*(1*(r8-R8)./tr8); 
    S19=S19+dt*(A(9)*(1-S19).*Hinf(V8-the(9))-B(9)*S19);
    
    %CNF alike PPN
    vcnf(:,i)=V9+dt*(1/Cm*(-Inal9-Ikl9-Ina9-Ik9-It9-Ihyp9-Inap9+Iappcnf-Igicnf-Isnrcnf));
    H9=H9+dt*((h9-H9)./th9);
    M9=M9+dt*((m9-M9)./tm9);
    Mk9=Mk9+dt*((mk9-Mk9)./tmk9);
    Mh9=Mh9+dt*((mh9-Mh9)./tmh9);
    Mt9=Mt9+dt*((mt9-Mt9)./tmt9);
    Ht9=Ht9+dt*((ht9-Ht9)./tht9);
    Hnap9=Hnap9+dt*((hnap9-Hnap9)./thnap9); 
    S18=S18+dt*(A(8)*(1-S18).*Hinf(V9-the(8))-B(8)*S18);
    
    %LC Arnauds paper
    vlc(:,i)=V10+dt*(1/Cm*(-Ina10-Ik10-Il10+Iapplc-Icnflc-Iprflc));
    q=q+dt*((qinf-q)./tq);
    S20=S20+dt*(A(10)*(1-S20).*Hinf(V10-the(10))-B(10)*S20);

end

%calculate freqeuncy for plotting
fr1=findfreq(vsn(1,:));
fr2=findfreq(vge(1,:));
fr3=findfreq(vgi(1,:));
fr4=findfreq(vppn(1,:));
fr5=findfreq(vth(1,:));
fr6=findfreq(vstr(1,:));
fr7=findfreq(vsnr(1,:));
fr8=findfreq(vprf(1,:));
fr9=findfreq(vcnf(1,:));
fr10=findfreq(vlc(1,:));

% [a,b,c]=periodogram(vsn(1,:),rectwin(length(vsn(1,:))),length(vsn(1,:)),100000,'ConfidenceLevel',0.95);
% figure
% plot(b,10*log10(c));
% xlim([13,30]);

%%Calculation of error index
GN=calculateEI(t,vth,timespike,tmax); %for thalamus
%GN=calculateEI(t,vppn,timespike,tmax); %for PPN

titleval=GN; %variable for plotting title

%%Plots membrane potential for one cell in each nucleus
%plotpotentials; 

return
