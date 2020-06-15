function DBS=createdbs(f,tmax,dt)
%Creates DBS train of frequency f, of length tmax (msec) with time step dt (msec)

t=0:dt:tmax; %create a time vector
DBS=zeros(1,length(t)); %create zeros for DBS
amp=300; %muAmp/cm^2
pulse=amp*ones(1,0.3/dt); %create a pulse with an amp and a width 0.3ms

i=1; %set a counter
%loop to create DBS train
while i<length(t)
    DBS(i:i+0.3/dt-1)=pulse; %create a pulse
    instfreq=f; %set frequency
    isi=tmax/instfreq; %inter-spike interval in ms
    i=i+round(isi/dt); %move i
end