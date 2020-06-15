function tau=ppn_taumt(V)
%tau=0.204+0.333./(exp(-(V+131)./16.7)+exp((V+15.8)./18.2));
tau=(0.612+1./(exp(-(V+134)./16.7)+exp((V+18.8)./18.2)))./6.9;
return