opt=load('optimization.mat');
opt2=load('optimization_str.mat');
h=opt.optimization(:,7);
pd=opt.optimization(:,8);
num=[];
for i=1:length(h)
    if h(i)<=0.04
        if pd(i)>=0.3
            num=[num,i];
        end
    end
end
        
        