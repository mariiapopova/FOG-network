function tau=snr_taumnap(V)
tau=0.03+0.116./(exp((-42.6-V)./14.4)+exp((-42.6-V)./-14.4));
return