function minf=snr_hinf(V)
minf=1./(1+exp((V+63.3)./8.1));
return