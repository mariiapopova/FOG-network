%function GN = FOGnetwork(pd,stim,freq, aopt, bopt, copt, dopt)
function GN = FOGnetwork(pd,stim,freq)

%Usage: GN = FOGnetwork(pd,stim,freq)
%
%Example: error_index=BG(1,1,130);
%
%Variables:
%pd - Variable to determine whether network is under the healthy or 
%Parkinsonian condition. For healthy, pd = 0, for Parkinson's, pd = 1.
%stim - Variable to determine whether deep brain stimulation is on.
%If DBS is off, stim = 1 If DBS is on, stim = 1.
%freq - Determines the frequency of stimulation, in Hz.
%
%
%Author: Mariia Popova, UKE; based on Rosa So, Duke University 
%Updated 15/06/2020

load('Istim.mat') %loads initial conditions %what is r?
addpath('gating')

%Membrane parameters
Cm=1; 

%Rm=100; Ie=0.7; Iege=0.2; Iegi=0; Iesnr=0.2; %wo syn
%Rm=100; Ie=0.26; Iege=1; Iegi=3; Iesnr=0.2; %wo syn
%Rm=100; Ie=1.5; Iege=5.6; Iegi=10; Iesnr=0.2; %wo syn -pd*2 
Rm=100; Ie=2; Iege=4.5; Iegi=8.5; Iesnr=1; %wo pd -pd*2 works!!! e2
Elsn=-70; Vth=-54; 
Vreset=-80; 

%In order of Th,STN,GP,PPN,Str or Th,STN,GPe,GPi,PPN,SNr,Str
gl=[0.05 2.25 0.1 0.1]; El=[-70 -60 -65 -67];
gna=[3 37 120 30 100]; Ena=[50 55 55 45 50]; 
gk=[5 45 30 3.2 80]; Ek=[-75 -80 -80 -95 -100];
gt=[5 0.5 0.5]; Et=0;
gnal5=0.0207; %na leakage ppn
gkl5=0.05; %k leakage ppn
gcort=0.15; Ecort=0; %cortex par for ppn
ghyp5=0.4; Ehyp5=-43; %hyperpolarization-activated current ppn
gnap5=1.9; %persistent sodium current ppn - chosen to have 8 Hz in rest
Pca=1*10^(-3); z=2; F=96490; R=8.314; T=309.5; Cao=2; Cai=2.4e-4; %alike destexhe ghk for ppn
Bcort=1;%ms^-1
gm=1; Em=-100;%for striatum muscarinic current

gca=[0 2 0.15]; Eca=[0 140 120];
gahp=[0 20 10]; %eahp and ek are the same excluding th
k1=[0 15 10]; %dissociation const of ahp current 
kca=[0 22.5 15]; %calcium pump rate constant

%synapse params alike in rubin SNr same to Gpi, PPN same to STN, Str same
%to Gpi
A=[0 5 2 2 5 2 2]; 
B=[0 1 0.04 0.08 1 0.08 0.08]; 
the=[0 30 20 20 30 20 20];

%%Synapse parameters
%In order of Igesn,Isnge,Igege,Isngi,Igegi,Igith 
Esyn = [-85 0 -85 0 -85 -85]; %alike in Rubin gsyn and in So Esyn
gsynppn = [0.26 0.2 0.26 0.2 0.2]; Esynppn = [0 -85 0 0 -85]; %in order snppn gippn ppnsn ppngi snrppn
tau=5; gpeak1=0.3; gpeak=0.43; %parameters for second-order alpha synapse
gsynsnr=0.3; Esynsnr=0; %for snr synapses in order stn
gsynstr=[0.5 0.5 0.5 0.5]; ggaba=0.1; gcorstr=0.07; Esynstr=[-85 0 -85 -85 -85 -85]; tau_i=13;%parameters for striatum synapses in order gaba-rec crtx strge gestr strgi strsnr

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

S5=zeros(n,1); %for alpha-synapse stn-ppn current
S6=zeros(n,1); %for alpha-synapse gpi-ppn current
S7=zeros(n,1); %for alpha-synapse ppn-stn current
S8=zeros(n,1); %for alpha-synapse ppn-gpi current
Sc=zeros(n,length(t)); %for cortex-ppn synapse

S1c=zeros(n,1); %for striatum gaba-rec
%for synapses striatum-gaba-rec
all=randsample(n,n);
bll=randsample(n,n);
cll=randsample(n,n);
dll=randsample(n,n);

S2a=zeros(n,1); 
S21a=zeros(n,1);
S2an=zeros(n,1);
S21an=zeros(n,1);
S2b=zeros(n,1); 
S21b=zeros(n,1);
S2c=zeros(n,1); 
S21c=zeros(n,1);
S3a=zeros(n,1); 
S31a=zeros(n,1); 
S3b=zeros(n,1);
S31b=zeros(n,1);
S32b=zeros(n,1);
S3c=zeros(n,1);
S31c=zeros(n,1); 
S32c=zeros(n,1);
S4=zeros(n,1); %for gpi-th current
S5=zeros(n,1);
S51=zeros(n,1);
S52=zeros(n,1);
S53=zeros(n,1);
S54=zeros(n,1);
S55=zeros(n,1);
S56=zeros(n,1);
S57=zeros(n,1);
S58=zeros(n,1);
S59=zeros(n,1);
S3d=zeros(n,1);
S9=zeros(n,1);
S91=zeros(n,1);
S92=zeros(n,1);
S93=zeros(n,1);
S94=zeros(n,1);
S95=zeros(n,1);
S96=zeros(n,1);
S97=zeros(n,1);
S98=zeros(n,1);
S99=zeros(n,1);

%%with or without DBS
Idbs=createdbs(freq,tmax,dt); %creating DBS train with frequency freq
if ~stim; Idbs=zeros(1,length(t)); end

%%initial conditions 
vth(:,1)=v1;
vsn(:,1)=v2;
vge(:,1)=v3;
vgi(:,1)=v4;
vppn(:,1)=v5; %for PPN
vsnr(:,1)=v6; %for SNr
vstr(:,1)=v7; %for striatum D2

%helper variables for gating and synapse params - starting parameters
N3=gpe_ninf(vge(:,1));%n for gpe
H3=gpe_hinf(vge(:,1));%h for gpe
R3=gpe_rinf(vge(:,1));%r for gpe
N4=gpe_ninf(vgi(:,1));%n for gpi
H4=gpe_hinf(vgi(:,1));%h for gpi
R4=gpe_rinf(vgi(:,1));%r for gpi
CA2=0.1; %intracellular concentration of Ca2+ in muM for stn
CA3=CA2; %for gpe
CA4=CA2; %for gpi
CA6=CA2; %for snr

H1=th_hinf(vth(:,1)); %h for th
R1=th_rinf(vth(:,1)); %r for th
M5=ppn_minf(vppn(:,1));%m for ppn na
H5=ppn_hinf(vppn(:,1));%h for ppn na
Mk5=ppn_mkinf(vppn(:,1));%m for ppn k
Mh5=ppn_mhinf(vppn(:,1));%m for ppn hyp
Hnap5=ppn_hnapinf(vppn(:,1));%h for ppn nap
Mt5=ppn_mtinf(vppn(:,1));%m for ppn low-tresh
Ht5=ppn_htinf(vppn(:,1));%h for ppn low-tresh

%striatum gating
m7=str_alpham(vstr(:,1))./(str_alpham(vstr(:,1))+str_betam(vstr(:,1)));
h7=str_alphah(vstr(:,1))./(str_alphah(vstr(:,1))+str_betah(vstr(:,1)));
n7=str_alphan(vstr(:,1))./(str_alphan(vstr(:,1))+str_betan(vstr(:,1)));
p7=str_alphap(vstr(:,1))./(str_alphap(vstr(:,1))+str_betap(vstr(:,1)));

timespikeint = int64(timespike); %index when 1 for ppn

%STN-GPe Synapse
const = gpeak/(tau*exp(-1));  
const1 = gpeak1/(tau*exp(-1)); 
const2 = gpeak1/(tau*exp(-1)); 
t_a = 1000; % Max duration of syn conductance
t_vec = 0:dt:t_a;

%STN-GPe Synapse
t_d_stn_gpe=2;
taudstngpea=2.5;
taurstngpea=0.4;
taudstngpen=67;
taurstngpen=2;
tpeakstngpea = t_d_stn_gpe + (((taudstngpea*taurstngpea)/(taudstngpea-taurstngpea))*log(taudstngpea/taurstngpea)); 
fstngpea = 1/(exp(-(tpeakstngpea-t_d_stn_gpe)/taudstngpea)-exp(-(tpeakstngpea-t_d_stn_gpe)/taurstngpea));
syn_func_stn_gpea = gpeak*fstngpea.*(exp(-(t_vec-t_d_stn_gpe)/taudstngpea)-exp(-(t_vec-t_d_stn_gpe)/taurstngpea)).*((t_vec>=t_d_stn_gpe)&(t_vec<=t_a));
tpeakstngpen = t_d_stn_gpe + (((taudstngpen*taurstngpen)/(taudstngpen-taurstngpen))*log(taudstngpen/taurstngpen)); 
fstngpen = 1/(exp(-(tpeakstngpen-t_d_stn_gpe)/taudstngpen)-exp(-(tpeakstngpen-t_d_stn_gpe)/taurstngpen));
syn_func_stn_gpen = gpeak*fstngpen.*(exp(-(t_vec-t_d_stn_gpe)/taudstngpen)-exp(-(t_vec-t_d_stn_gpe)/taurstngpen)).*((t_vec>=t_d_stn_gpe)&(t_vec<=t_a));

%STN-GPi Synapse %SNr in a similar fashion
t_d_stn_gpi=1.5;
syn_func_stn_gpi = const*(t_vec-t_d_stn_gpi).*(exp(-(t_vec-t_d_stn_gpi)/tau)).*((t_vec>=t_d_stn_gpi)&(t_vec<=t_a));

%GPe-STN Synapse
t_d_gpe_stn=4;
taudg=7.7;
taurg=0.4;
tpeakg = t_d_gpe_stn + (((taudg*taurg)/(taudg-taurg))*log(taudg/taurg)); 
fg = 1/(exp(-(tpeakg-t_d_gpe_stn)/taudg)-exp(-(tpeakg-t_d_gpe_stn)/taurg));
syn_func_gpe_stn = gpeak1*fg.*(exp(-(t_vec-t_d_gpe_stn)/taudg)-exp(-(t_vec-t_d_gpe_stn)/taurg)).*((t_vec>=t_d_gpe_stn)&(t_vec<=t_a));

%GPe-GPi Synapse
t_d_gpe_gpi=3;
syn_func_gpe_gpi = const1*(t_vec-t_d_gpe_gpi).*(exp(-(t_vec-t_d_gpe_gpi)/tau)).*((t_vec>=t_d_gpe_gpi)&(t_vec<=t_a));

%GPe-GPe Synapse
t_d_gpe_gpe=1;
syn_func_gpe_gpe = const1*(t_vec-t_d_gpe_gpe).*(exp(-(t_vec-t_d_gpe_gpe)/tau)).*((t_vec>=t_d_gpe_gpe)&(t_vec<=t_a));

%GPi-TH Synapse
t_d_gpi_th=5;
syn_func_gpi_th = const1*(t_vec-t_d_gpi_th).*(exp(-(t_vec-t_d_gpi_th)/tau)).*((t_vec>=t_d_gpi_th)&(t_vec<=t_a));

%Str-GPe Synapse
t_d_d2_gpe=5;
syn_func_str_gpe = const2*(t_vec- t_d_d2_gpe).*(exp(-(t_vec-t_d_d2_gpe)/tau)).*((t_vec>=t_d_d2_gpe)&(t_vec<=t_a));

%Str-GPi Synapse
t_d_d1_gpi=4;
syn_func_str_gpi = const2*(t_vec- t_d_d1_gpi).*(exp(-(t_vec-t_d_d1_gpi)/tau)).*((t_vec>=t_d_d1_gpi)&(t_vec<=t_a));

t_list_str_indr(1:n) = struct('times',[]);
t_list_stn(1:n) = struct('times',[]);
t_list_gpe(1:n) = struct('times',[]);
t_list_gpi(1:n) = struct('times',[]);

ggege=0.5;
gcorsna=0.3;
gcorsnn=0.003;
gsngen=0.001;
gsngea=0.15;
gsngi=0.15;
ggith=0.112;
ggesn=0.5;
gstrgpe=0.5;
gstrgpi=0.5;
ggigi=0.5;

%Time loop
for i=1:length(t)-1
    
    %condition for cortex current for ppn and striatum
    if ismember(i+1, timespikeint/dt)
        Sc(:,i) = 1;
    end
    
    V1=vth(:,i);    V2=vsn(:,i);     V3=vge(:,i);    V4=vgi(:,i); V5=vppn(:,i); V6=vsnr(:,i); 
    V7=vstr(:,i);

    % Synapse parameters 
    S21a(2:n)=S2a(1:n-1);
    S21a(1)=S2a(n);
    
    S21an(2:n)=S2an(1:n-1);
    S21an(1)=S2an(n);
    
    S21b(2:n)=S2b(1:n-1);
    S21b(1)=S2b(n);
    
    S21c(2:n)=S2c(1:n-1);
    S21c(1)=S2c(n);
    
    S31a(1:n-1)=S3a(2:n);
    S31a(n)=S3a(1);
    
    S31b(1:n-1)=S3b(2:n);
    S31b(n)=S3b(1);
    
    S31c(1:n-1)=S3c(2:n);
    S31c(n)=S3c(1);
    
    S32c(3:n)=S3c(1:n-2);
    S32c(1:2)=S3c(n-1:n);
    
    S32b(3:n)=S3b(1:n-2);
    S32b(1:2)=S3b(n-1:n);
    
    S51(1:n-1)=S5(2:n);
    S51(n)=S5(1);
    S52(1:n-2)=S5(3:n);
    S52(n-1:n)=S5(1:2);
    S53(1:n-3)=S5(4:n);
    S53(n-2:n)=S5(1:3);
    S54(1:n-4)=S5(5:n);
    S54(n-3:n)=S5(1:4);
    S55(1:n-5)=S5(6:n);
    S55(n-4:n)=S5(1:5);
    S56(1:n-6)=S5(7:n);
    S56(n-5:n)=S5(1:6);
    S57(1:n-7)=S5(8:n);
    S57(n-6:n)=S5(1:7);
    S58(1:n-8)=S5(9:n);
    S58(n-7:n)=S5(1:8);
    S59(1:n-9)=S5(10:n);
    S59(n-8:n)=S5(1:9);
    
    S91(1:n-1)=S9(2:n);
    S91(n)=S9(1);
    S92(1:n-2)=S9(3:n);
    S92(n-1:n)=S9(1:2);
    S93(1:n-3)=S9(4:n);
    S93(n-2:n)=S9(1:3);
    S94(1:n-4)=S9(5:n);
    S94(n-3:n)=S9(1:4);
    S95(1:n-5)=S9(6:n);
    S95(n-4:n)=S9(1:5);
    S96(1:n-6)=S9(7:n);
    S96(n-5:n)=S9(1:6);
    S97(1:n-7)=S9(8:n);
    S97(n-6:n)=S9(1:7);
    S98(1:n-8)=S9(9:n);
    S98(n-7:n)=S9(1:8);
    S99(1:n-9)=S9(10:n);
    S99(n-8:n)=S9(1:9);
    
    S11cr=S1c(all); S12cr=S1c(bll); S13cr=S1c(cll); S14cr=S1c(dll); %dummy striatum gaba current
    
    %membrane parameters - gating variables
    m1=th_minf(V1);m3=gpe_minf(V3);m4=gpe_minf(V4); m6=gpe_minf(V6); %gpe and gpi are modeled similarily
    n3=gpe_ninf(V3);n4=gpe_ninf(V4); n6=gpe_ninf(V6);
    h1=th_hinf(V1);h3=gpe_hinf(V3);h4=gpe_hinf(V4); h6=gpe_hinf(V6);
    p1=th_pinf(V1);
    a3=gpe_ainf(V3);a4=gpe_ainf(V4); a6=gpe_ainf(V6); %for low-treshold ca
    s3=gpe_sinf(V3);s4=gpe_sinf(V4); s6=gpe_sinf(V6);
    r1=th_rinf(V1);r3=gpe_rinf(V3);r4=gpe_rinf(V4); r6=gpe_rinf(V6);
    m5=ppn_minf(V5); h5=ppn_hinf(V5); mk5=ppn_mkinf(V5); mh5=ppn_mhinf(V5);
    mnap5=ppn_mnapinf(V5); hnap5=ppn_hnapinf(V5); mt5=ppn_mtinf(V5); ht5=ppn_htinf(V5);

    %membrane parameters - time constants
    tn3=gpe_taun(V3);tn4=gpe_taun(V4);tn6=gpe_taun(V6);
    th1=th_tauh(V1);th3=gpe_tauh(V3);th4=gpe_tauh(V4); th6=gpe_tauh(V6);
    tr1=th_taur(V1);tr3=30;tr4=30;tr6=30;
    tm5=ppn_taum(V5); th5=ppn_tauh(V5); tmk5=ppn_taumk(V5); tmh5=ppn_taumh(V5);
    thnap5=ppn_tauhnap(V5); tmt5=ppn_taumt(V5); tht5=ppn_tauht(V5);
    
    %thalamic cell currents
    Il1=gl(1)*(V1-El(1));
    Ina1=gna(1)*(m1.^3).*H1.*(V1-Ena(1));
    Ik1=gk(1)*((0.75*(1-H1)).^4).*(V1-Ek(1)); %misspelled in So paper
    It1=gt(1)*(p1.^2).*R1.*(V1-Et);
    Igith=ggith*(V1-Esyn(6)).*(S4); 
    
    %STN cell currents
    Igesn=(ggesn*(V2-Esyn(1)).*(S3a+S31a)); 
    Icorsnampa=gcorsna.*(V2-Esynstr(2)).*(2*Sc(:,i));
    Icorsnnmda=gcorsnn.*(V2-Esynstr(2)).*(2*Sc(:,i));
    Ippnsn=gsynppn(3)*(V2-Esynppn(3)).*S7; %first-order kinetics ppn to stn
    
    %GPe cell currents
    Il3=gl(3)*(V3-El(3));
    Ik3=gk(3)*(N3.^4).*(V3-Ek(3));
    Ina3=gna(3)*(m3.^3).*H3.*(V3-Ena(3));
    It3=gt(3)*(a3.^3).*R3.*(V3-Eca(3)); %Eca as in Rubin and Terman
    Ica3=gca(3)*(s3.^2).*(V3-Eca(3)); %misspelled in So paper
    Iahp3=gahp(3)*(V3-Ek(3)).*(CA3./(CA3+k1(3))); %as Ek is the same with Eahp
    Isnge=(gsngea)*((V3-Esyn(2)).*(S2a+S21a)); 
    Isngenmda=(gsngen)*((V3-Esyn(2)).*(S2an+S21an)); 
    %Igege=0.25*((ggege+3*pd*ggege).*(V3-Esyn(3)).*(S31c+S32c));
    %Igege=0.25*((ggege+5*pd*ggege).*(V3-Esyn(3)).*(S31c+S32c)); smh
    %Igege=0.25*((ggege+dopt*pd*ggege).*(V3-Esyn(3)).*(S31c+S32c));
    Igege=0.25*((ggege+9*pd*ggege).*(V3-Esyn(3)).*(S31c+S32c));
    %str-gpe synapse
    Istrge=gstrgpe*(V3-Esynstr(3)).*(S5+S51+S52+S53+S54+S55+S56+S57+S58+S59);

    %GPi cell currents
    Il4=gl(3)*(V5-El(3));
    Ik4=gk(3)*(N4.^4).*(V4-Ek(3));
    Ina4=gna(3)*(m4.^3).*H4.*(V4-Ena(3)); %Eca as in Rubin and Terman
    It4=gt(3)*(a4.^3).*R4.*(V4-Eca(3)); %misspelled in So paper
    Ica4=gca(3)*(s4.^2).*(V4-Eca(3)); 
    Iahp4=gahp(3)*(V4-Ek(3)).*(CA4./(CA4+k1(3))); %as Ek is the same with Eahp
    Isngi=gsngi*(V4-Esyn(4)).*(S2b+S21b); 
    Igegi=ggigi*(V4-Esyn(5)).*(S31b+S32b); 
    Ippngi=gsynppn(4)*(V4-Esynppn(4)).*S8; %first-order kinetics ppn to gpi
    %str-gpi synapse
    Istrgi=gstrgpi*(V4-Esynstr(5)).*(S9+S91+S92+S93+S94+S95+S96+S97+S98+S99);
    
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
    Icort=gcort*Sc(:,i).*(V5-Ecort);
    %Isnrppn=gsynppn(5)*(V5-Esynppn(5)).*S14; %first-order kinetics snr to ppn
    
    %SNr cell currents - modelled as GPi
    Isnsnr=gsngi.*(V6-Esynsnr).*(S2c+S21c); 
    %str-snr synapse
    %Istrsnr=gsynstr(4)*(V6-Esynstr(6)).*S13; %1str to 1ge
    
    %Striatum D2 cell currents
    Ina7=gna(5)*(m7.^3).*h7.*(V7-Ena(5));
    Ik7=gk(5)*(n7.^4).*(V7-Ek(5));
    Il7=gl(4)*(V7-El(4));
    %Im7=(2.6-1.1*pd)*gm*p7.*(V7-Em);
    %Im7=(2.6-1.1*pd)*gm*p7.*(V7-Em); smh
    %Im7=(2.6-copt*pd)*gm*p7.*(V7-Em);
    Im7=(2.6-0.1*pd)*gm*p7.*(V7-Em);
    Igaba7=(ggaba/4)*(V7-Esynstr(1)).*(S11cr+S12cr+S13cr+S14cr); %maybe change to 3.5 for direct and indirect %recieves input from 40% remaining
    %Icorstr= (gcorstr-0.044*pd)*(V7-Esynstr(2)).*Sc(:,i);
    %Icorstr=(5*gcorstr-0.044*pd)*(V7-Esynstr(2)).*Sc(:,i); smh
    %Icorstr=(aopt*gcorstr-bopt*pd)*(V7-Esynstr(2)).*Sc(:,i);
    Icorstr=(4.75*gcorstr-0.154*pd)*(V7-Esynstr(2)).*Sc(:,i);
    %ge-str synapse
    Igestr=gstrgpe*(V7-Esynstr(3)).*(S3d);
    Iappstr=0; 
    
    %Differential Equations for cells using forward Euler method
    %thalamic
    vth(:,i+1)= V1+dt*(1/Cm*(-Il1-Ik1-Ina1-It1-Igith+Istim(i+1)));
    H1=H1+dt*((h1-H1)./th1);
    R1=R1+dt*((r1-R1)./tr1);
    Sc(:,i+1)=Sc(:,i)+dt*(-Bcort.*Sc(:,i));
    
    %STN 
    vsn(:,i+1) = V2 + (dt/Cm)*(((Elsn-V2)/Rm) + Ie -Icorsnampa-Icorsnnmda-Igesn);%+Idbs(i)-Ippnsn); %currently STN-DBS
    
    if vsn(:,i+1)>= Vth
        vsn(:,i)=40;
        vsn(:,i+1)=Vreset;
    end 
    
    
    %Calculate synapses
    for j=1:n
        if (vsn(j,i+1)<-10 && vsn(j,i)>-10) % check for input spike
             t_list_stn(j).times = [t_list_stn(j).times; 1];
        end   
        
       %Calculate synaptic current due to current and past input spikes
       S2a(j) = sum(syn_func_stn_gpea(t_list_stn(j).times));
       S2an(j) = sum(syn_func_stn_gpen(t_list_stn(j).times));
       S2b(j) = sum(syn_func_stn_gpi(t_list_stn(j).times));
       S2c(j) = sum(syn_func_stn_gpi(t_list_stn(j).times));

       %Update spike times
       if t_list_stn(j).times
         t_list_stn(j).times = t_list_stn(j).times + 1;
         if (t_list_stn(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
           t_list_stn(j).times = t_list_stn(j).times((2:max(size(t_list_stn(j).times))));
         end
       end
     end

    %GPe
%     vge(:,i+1) = V3 + (dt/Cm)*(((Elsne-V3)/Rm) + Iege-Isnge-Igege-Isngenmda-Istrge);
%     
%     if vge(:,i+1)>= Vth
%         V3=40;
%         vge(:,i+1)=Vresete;
%     end 
%     vge(:,i)=V3;
    
    vge(:,i+1)=V3+dt*(1/Cm*(-Il3-Ik3-Ina3-It3-Ica3-Iahp3-Isnge-Isngenmda-Igege+Iege-Istrge));
    
    N3=N3+dt*(0.1*(n3-N3)./tn3); %misspelled in So paper
    H3=H3+dt*(0.05*(h3-H3)./th3); %misspelled in So paper
    R3=R3+dt*(1*(r3-R3)./tr3); %misspelled in So paper
    CA3=CA3+dt*(1*10^-4*(-Ica3-It3-kca(3)*CA3));
    
    %Calculate synapses
    for j=1:n
        if (vge(j,i)<-10 && vge(j,i+1)>-10) % check for input spike
             t_list_gpe(j).times = [t_list_gpe(j).times; 1];
        end  
        
       %Calculate synaptic current due to current and past input spikes
       S3a(j) = sum(syn_func_gpe_stn(t_list_gpe(j).times));
       S3b(j) = sum(syn_func_gpe_gpi(t_list_gpe(j).times));
       S3c(j) = sum(syn_func_gpe_gpe(t_list_gpe(j).times));
       S3d(j) = sum(syn_func_str_gpe(t_list_gpe(j).times));
   
       %Update spike times
       if t_list_gpe(j).times
         t_list_gpe(j).times = t_list_gpe(j).times + 1;
         if (t_list_gpe(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
           t_list_gpe(j).times = t_list_gpe(j).times((2:max(size(t_list_gpe(j).times))));
         end
       end
    end
    
    %GPi
%     vgi(:,i+1) = V4 + (dt/Cm)*(((Elsne-V4)/Rm) + Iegi-Isngi-Igegi-Istrgi);%-Ippngi;
%     
%     if vgi(:,i+1)>= Vth
%         V4=40;
%         vgi(:,i+1)=Vreset;
%     end 
%     vgi(:,i)=V4;
    
    vgi(:,i+1)=V4+dt*(1/Cm*(-Il4-Ik4-Ina4-It4-Ica4-Iahp4+Iegi-Isngi-Igegi-Istrgi));%-Ippngi
    
    N4=N4+dt*(0.1*(n4-N4)./tn4); %misspelled in So paper
    H4=H4+dt*(0.05*(h4-H4)./th4); %misspelled in So paper
    R4=R4+dt*(1*(r4-R4)./tr4); %misspelled in So paper
    CA4=CA4+dt*(1*10^-4*(-Ica4-It4-kca(3)*CA4));
    
    %Calculate synapses
    for j=1:n
        if (vgi(j,i)<-10 && vgi(j,i+1)>-10) % check for input spike
             t_list_gpi(j).times = [t_list_gpi(j).times; 1];
        end
        
        % Calculate synaptic current due to current and past input spikes
        S4(j) = sum(syn_func_gpi_th(t_list_gpi(j).times));


       % Update spike times
       if t_list_gpi(j).times
         t_list_gpi(j).times = t_list_gpi(j).times + 1;
         if (t_list_gpi(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
           t_list_gpi(j).times = t_list_gpi(j).times((2:max(size(t_list_gpi(j).times))));
         end
       end
    end
    
    %PPN
    vppn(:,i+1)=V5+dt*(1/Cm*(-Inal5-Ikl5-Ina5-Ik5-It5-Ihyp5-Inap5+Iappppn-Icort));%-Isnppn-Igippn-Isnrppn
    
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
    
    %SNr
    %vsnr(:,i)=V6+dt*(1/Cm*(-Il6-Ik6-Ina6-It6-Ica6-Iahp6-Isnsnr-Istrsnr+Idbs(i)+Iappsnr)); 
    vsnr(:,i+1)=V6 + (dt/Cm)*(((Elsn-V6)/Rm) + Iesnr);%-Isnsnr-Istrsnr
    
    if vsnr(:,i+1)>= Vth
        V6=40;
        vsnr(:,i+1)=Vreset;
    end 
    vsnr(:,i)=V6;
     
    %S14=S14+dt*(A(6)*(1-S14).*Hinf(V6-the(6))-B(6)*S14); %snr-ppn synapse

    
    %Striatum D2
    vstr(:,i+1)=V7+(dt/Cm)*(-Ina7-Ik7-Il7-Im7-Igaba7-Icorstr+Iappstr-Igestr);
    
    m7=m7+dt*(str_alpham(V7).*(1-m7)-str_betam(V7).*m7);
    h7=h7+dt*(str_alphah(V7).*(1-h7)-str_betah(V7).*h7);
    n7=n7+dt*(str_alphan(V7).*(1-n7)-str_betan(V7).*n7);
    p7=p7+dt*(str_alphap(V7).*(1-p7)-str_betap(V7).*p7);
    S1c=S1c+dt*((str_Ggaba(V7).*(1-S1c))-(S1c/tau_i));
    
    for j=1:n
        if (vstr(j,i)<-10 && vstr(j,i+1)>-10) % check for input spike
             t_list_str_indr(j).times = [t_list_str_indr(j).times; 1];
        end  
        
       % Calculate synaptic current due to current and past input spikes
       S5(j) = sum(syn_func_str_gpe(t_list_str_indr(j).times));
       S9(j) = sum(syn_func_str_gpi(t_list_str_indr(j).times));

       % Update spike times
       if t_list_str_indr(j).times
         t_list_str_indr(j).times = t_list_str_indr(j).times + 1;
         if (t_list_str_indr(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
            t_list_str_indr(j).times = t_list_str_indr(j).times((2:max(size(t_list_str_indr(j).times))));
         end
       end
    end
    
end

%calculate freqeuncy for plotting
fr1=findfreq(vsn(1,:));
fr2=findfreq(vge(1,:));
fr3=findfreq(vgi(1,:));
fr4=findfreq(vppn(1,:));
fr5=findfreq(vth(1,:));
fr6=findfreq(vstr(1,:));
fr7=findfreq(vsnr(1,:));

%%Calculation of error index
GN=calculateEI(t,vth,timespike,tmax); %for thalamus
%GN=calculateEI(t,vppn,timespike,tmax); %for PPN

titleval=GN; %variable for plotting title

%Plots membrane potential for one cell in each nucleus
plotpotentials; 

return
