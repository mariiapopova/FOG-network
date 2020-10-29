function y = SERmodel_multneuro(NN,C,T,ia,dbs,d)

% SER network simulation for multiple neurons in one region
% 
% NN           = Number of neurons in a region
% C            = Matrix of coupling (NxN) between pairs of regions (can be directed) 
% T            = Total time of simulated activity
% ia           = Initial condition setting. It can be a number representing the number of excited nodes (the remaining nodes are splitted in two equal size cohorts of susceptible and refractory nodes)
%                                           or can be a vector describing
%                                           the initial state of each
%                                           region
% d            = Initial time steps to remove (transient dynamics)
% 
% Convention is:
%       - susceptible node =  0
%       - excited node     =  1
%       - refractory node  = -1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6, d = 0; end
if nargin<5, ia = 1; end
if T<=d, error(strcat('Simulation time must be greater than the transient ',num2str(d))); end
if ia>length(C), error('Initial active nodes must be equal to or lower than the total number of nodes'); end

numb=NN; %number of neurons in region
N   = length(C);
y   = zeros(numb,N,T,'single');     % initialize phase timeseries for one cycle

% Initialization

for j = 1:numb
    if length(ia)==N
        y(j,:,1) = ia;
    else
        disp('here');
        y(j,randsample(N,ia),1) = 1;
        id = find(y(j,:,1)==0);
        y(j,id,1) = floor(1-2*rand(length(id),1));
    end
end

% Equations integration

for t = 1:T-1
    for i = 1:numb
        y(i, y(i,:,t)==1,t+1) = 1;
        y(i, y(i,:,t)==1 & (rand(1,10)>0.5),t+1) = -1;
        y(i,10,t+1) = 1; %cortex always
        if dbs==1
            y(i,6,t+1) = 1; %stn always
        end
        if dbs==2
            y(i,6,t+1) = 1; %stn always
            y(i,5,t+1) = 1; %snr always
        end
        %y(i, y(i,:,t)==1 & (sum(C(:,y(i,:,t)==1),2)>0).' & (rand(1,10)>0.5),t+1) = 1;
        y(i, y(i,:,t)==0 & (sum(C(:,y(i,:,t)==1),2)>0).' & (rand(1,10)>0.5),t+1) = 1;
    end
end

y(:,:,1:d) = [];  % remove initial dis steps of simulations to exclude transient dynamics