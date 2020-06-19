function minf=snr_snainf(V)
minf=1./(1+exp((V+30)./0.4));
return