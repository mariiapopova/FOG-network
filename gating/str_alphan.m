function an = str_alphan(V)
an=(0.032*(52+V))./(1-exp((-52-V)./5));
end