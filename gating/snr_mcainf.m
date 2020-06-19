function minf=snr_mcainf(V)
minf=1./(1+exp(-(V+27.5)./3));
return