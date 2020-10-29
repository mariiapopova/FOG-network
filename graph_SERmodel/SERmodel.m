function y = Network_SER(C,T,ia,d)

% SER network simulation for 1 model in a region
% 
% C            = Matrix of coupling (NxN) between pairs of regions (can be directed) 
% T            = Total time of simulated activity
% ia           = Initial condition setting. It can be a number representing the number of excited nodes (the remaining nodes are splitted in two equal size cohorts of susceptible and refractory nodes)
%                                           or can be a vector describing the initial state of each node
% d            = Initial time steps to remove (transient dynamics)
% 
% Convention is:
%       - susceptible node =  0
%       - excited node     =  1
%       - refractory node  = -1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4, d = 0; end
if nargin<3, ia = 1; end
if T<=d, error(strcat('Simulation time must be greater than the transient ',num2str(d))); end
if ia>length(C), error('Initial active nodes must be equal to or lower than the total number of nodes'); end

N   = length(C);
y   = zeros(N,T,'single');     % initialize phase timeseries for one cycle

% Initialization

if length(ia)==N
    y(:,1) = ia;
else
    y(randsample(N,ia),1) = 1;
    id = find(y(:,1)==0);
    y(id,1) = floor(1-2*rand(length(id),1));
end

% Equations integration

for t = 1:T-1
    y(y(:,t)==1,t+1) = -1;
    y(10,t+1) = 1; %cortex always
    %y(6,t) = 1; %stn always
    %y(5,t) = 1; %snr always
    %y(y(:,t)==1 & sum(C(:,y(:,t)==1),2)>0,t+1) = 1;
    y(y(:,t)==0 & sum(C(:,y(:,t)==1),2)>0,t+1) = 1;
end

y(:,1:d) = [];  % remove initial dis steps of simulations to exclude transient dynamics
