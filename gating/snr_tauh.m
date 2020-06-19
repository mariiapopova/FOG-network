function tau=snr_tauh(V)
tau=0.59+34.5./(exp((-43-V)./10)+exp((-43-V)./-5));
return