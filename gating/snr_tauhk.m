function tau=snr_tauhk(V)
tau=5+15./(exp((-V)./10)+exp((-V)./-10));
return