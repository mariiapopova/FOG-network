function minf=lc_minf(V)
a=0.1*(V+29.7)./(1-exp(-(V+29.7)./10));
b=4.*exp(-(V+54.7)./18);
minf=a./(a+b);
return