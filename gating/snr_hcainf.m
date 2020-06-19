function minf=snr_hcainf(V)
minf=1./(1+exp((V+52.5)./5.2));
return