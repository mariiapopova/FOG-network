function minf=snr_minf(V)
minf=1./(1+exp(-(V+30.2)./6.2));
return