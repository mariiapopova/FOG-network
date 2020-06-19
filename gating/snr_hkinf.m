function minf=snr_hkinf(V)
minf=1./(1+exp((V+20)./10));
return