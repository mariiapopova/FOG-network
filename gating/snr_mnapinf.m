function minf=snr_mnapinf(V)
minf=1./(1+exp(-(V+50)./3));
return