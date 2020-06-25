clear all
close all
addpath('gating')
%check why the form is a lil offf!

Cm=1; 

gna=[3 37 120 30]; Ena=[50 55 55 45]; 
gk=[5 45 30 3.2]; Ek=[-75 -80 -80 -95];
gnal5=0.0207; %na leakage ppn
gkl5=0.05; %k leakage ppn
ghyp5=0.4; Ehyp5=-43; %hyperpolarization-activated current ppn
gnap5=1.9; %persistent sodium current ppn
Pca=1*10^(-3); z=2; F=96490; R=8.314; T=309.5; Cao=2; Cai=2.4e-4; %alike destexhe ghk for ppn

%time variables
tmax=1000; %maximum time (ms)
dt=0.01; %timestep (ms)
t=0:dt:tmax; %time vector
n=10; %number of neurons in each nucleus (TH, STN, GPe, GPi)

v5=-62+randn(n,1)*5; %for PPN
vppn=zeros(n,length(t)); %PPN membrane voltage
vppn(:,1)=v5; %for PPN

M5=ppn_minf(vppn(:,1));%m for ppn na
H5=ppn_hinf(vppn(:,1));%h for ppn na
Mk5=ppn_mkinf(vppn(:,1));%m for ppn k
Mh5=ppn_mhinf(vppn(:,1));%m for ppn hyp
Hnap5=ppn_hnapinf(vppn(:,1));%h for ppn nap
Mt5=ppn_mtinf(vppn(:,1));%m for ppn low-tresh
Ht5=ppn_htinf(vppn(:,1));%h for ppn low-tresh

for i=2:length(t)  
    V5=vppn(:,i-1);
    m5=ppn_minf(V5); h5=ppn_hinf(V5); mk5=ppn_mkinf(V5); mh5=ppn_mhinf(V5);
    mnap5=ppn_mnapinf(V5); hnap5=ppn_hnapinf(V5); mt5=ppn_mtinf(V5); ht5=ppn_htinf(V5);
    tm5=ppn_taum(V5); th5=ppn_tauh(V5); tmk5=ppn_taumk(V5); tmh5=ppn_taumh(V5);
    thnap5=ppn_tauhnap(V5); tmt5=ppn_taumt(V5); tht5=ppn_tauht(V5);

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
    
    %PPN
    vppn(:,i)=V5+dt*(1/Cm*(-Inal5-Ikl5-Ina5-Ik5-It5-Ihyp5-Inap5+Iappppn));
    H5=H5+dt*((h5-H5)./th5);
    M5=M5+dt*((m5-M5)./tm5);
    Mk5=Mk5+dt*((mk5-Mk5)./tmk5);
    Mh5=Mh5+dt*((mh5-Mh5)./tmh5);
    Mt5=Mt5+dt*((mt5-Mt5)./tmt5);
    Ht5=Ht5+dt*((ht5-Ht5)./tht5);
    Hnap5=Hnap5+dt*((hnap5-Hnap5)./thnap5); %alike rubin
end

figure
plot(vppn(1,:));
%% kostya
% figure
% dt=1e-5;
% val = vppn(1,:);
% val_mean0 = val - mean(val)
% Fs = 1/dt;                   % samples per second
% N = length(val);            % samples
% dF = Fs/N;                 % hertz per sample
% f = -Fs/2:dF:Fs/2-dF + (dF/2)*mod(N,2);      % hertz
% val_fft_centered = fftshift(fft(val_mean0))/N;
% figure;
% plot(f,abs(val_fft_centered));

fr=findfreq(vppn(1,:));
disp(fr);

