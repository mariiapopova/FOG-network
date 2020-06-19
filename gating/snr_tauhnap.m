function tau=snr_tauhnap(V)
tau=10+7./(exp((-34-V)./26)+exp((-34-V)./-31.9));
return