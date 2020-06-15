function [Istim, timespike]= createSMCinput(tmax,dt,freq,cv)
%creates Sensorimotor Cortex (SMC) input to thalamic cells

%Variables:
%tmax - length of input train (msec)
%dt - time step (msec)
%freq - frequency of input train
%cv - coefficient of variation of input train (gamma distribution)

%Output
%Istim - Input train from SMC
%timespike - Timing of each input pulse

t=0:dt:tmax; %create time vector
amp=3.5; %pulse amplitude in muA/cm^2
Istim=zeros(1,length(t)); %empty stim array
dur=5; %pulse duration in ms
pulse=amp*ones(1,dur/dt); %create a pulse

i=1; j=1; %set counters
A = 1/cv^2; %shape parameter 
B  = freq / A; %scale parameter 
if cv==0
    instfreq=freq; %take freq
else
    instfreq=gamrnd(A,B); %take freq from gamma distribution
end

ipi=tmax/instfreq; %calculate inter-spike interval!!! was 1000 before
i=i+round(ipi/dt); %inter-spike interval

%loop to create SMC input
while i<length(t) 
    timespike(j)=t(i); %calculate time when spike happens (begins)
    Istim(i:i+dur/dt-1)=pulse; %create stimulus spike
    %A = 1/cv^2;
    %B  = freq / A;
    %to create variability in ipi
    if cv==0
        instfreq=freq;
    else
        instfreq=gamrnd(A,B);
    end
    ipi=tmax/instfreq;
    i=i+round(ipi/dt); %move for one ipi
    j=j+1; %increase j
end

return



