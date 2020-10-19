function binf=lc_binf(V)
binf=(1./(1+exp(0.069.*(V+53.3)))).^4;
return