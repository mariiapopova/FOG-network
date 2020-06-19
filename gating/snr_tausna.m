function tau=snr_tausna(V)
tau=10+40./(exp((-40-V)./18.3)+exp((-40-V)./-10));
return