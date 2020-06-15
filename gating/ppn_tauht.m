function tau=ppn_tauht(V)
if (V+2)<-81
	%tau=0.333.*exp((V+466)./66.6);
    tau=exp((V+469)./66.6)./3.74;
else
    %tau=9.32+0.333.*exp(-(V+21)./10.5);
    tau=(28+exp(-(V+24)./10.5))./3.74;
end
return