opt=load('optimization.mat');
opt2=load('optimization_str.mat');
h=opt.optimization(:,4);
pd=opt.optimization(:,5);
num=[];
for i=1:length(h)
    %if h(i)<=0.01
        if pd(i)>=0.5
            num=[num,i];
        end
    %end
end
        
        