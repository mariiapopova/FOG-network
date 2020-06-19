function tau=snr_taumk(V)
tau=0.1+13.9./(exp((-26-V)./13)+exp((-26-V)./-12));
return