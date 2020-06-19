function minf=snr_mkinf(V)
minf=1./(1+exp(-(V+26)./7.8));
return