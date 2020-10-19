function qinf=lc_qinf(V)
a=0.01.*(V+45.7)./(1-exp(-(V+45.7)./10));
b=0.125.*exp(-(V+55.7)./80);
ninf=a./(a+b);
binf=(1./(1+exp(0.069.*(V+53.3)))).^4;
qinf=ninf.^4+(0.21*(47.7/20)).*binf;
return