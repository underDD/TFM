% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons

%% SET UP THE NETWORK

close all; clear all;

% 1- NETWORK SIZE:
Ne=800; Ni=200; % Excitatory, inhibitory. Ne+Ni is total neurons.
SIM_TIME=3000;

% 2 - GLOBAL PARAMETERS THAT SET OUR NEURON MODEL. DEFAULT IS SPIKING
% NEURON:
% Set initial conditions of neurons, with some variability provided by the
% vectors re and ri containing random numbers between 0 and 1.
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri]; % Model parameters. Do not touch.
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];

% 3 - SET UP THE CONNECTIVITY MATRIX: DIRECTED NETWORK
% In this construction, 1=connection exist, 0=no connection.
% Connectivity is set as random. Then, a fraction of connections are set 0.
% Note that effectively this is an Erdös-Rényi graphs, with no spatial
% characteristics. % If you want to symmetrize the network to make 
% undirected networks (not realistic in neuroscience), use A = (A + A.')/2
% after line 27.
frac_delete=0.8; % Set this fraction of connections to zero.
A=[rand(Ne+Ni)]; % Set a random graph by filling the matrix with random 
                 % numbers [0,1].
A(A<frac_delete)=0; % Eliminate connections to have a sparse network.
A(A>0)=1; % Set surviving connections to 1 (binary matrix).     
A = A - diag(diag(A)); % Make the diagonal elements 0 (no self-connection).


figure;
imagesc(A);
xlim([1,1000]); xlabel('neuron');
ylim([1,1000]); ylabel('neuron');
title('Connectivity matrix');
axis tight;
axis equal;

% 4 - SET SYNAPTIC WEIGHTS (STRENGTHS) OF CONNECTIONS.
% I represents the amnplitudes of evoked currents (EPSC and IPSC). 
% Increase the max. excitatory weigth to induce synchronization.
% Be careful. The balance between excitation and inhibition governs the
% dynamics of the network!
MAX_EXC_WEIGTH=4; MAX_INH_WEIGTH=.5; 
W=[MAX_EXC_WEIGTH*rand(Ne+Ni,Ne), -MAX_INH_WEIGTH*rand(Ne+Ni,Ni)];

% 5 - The final connectivity matrix S is the element-wise multiplication of A
% and W.
S=A.*W;  % Note that S is directed and weighted!!

% 6 - DEFINE NOISE STRENGTH. Remember that noise (or random inputs) is the
% main drive of spontaneous activity. Default is around 5. 
NOISE_MAX=3;


%% MAIN SIMULATION:

v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u
firings=[]; % spike timings
for t=1:SIM_TIME % simulation of 1000 ms
I=[NOISE_MAX*randn(Ne,1);2*randn(Ni,1)]; % NOISE or thalamic input
fired=find(v>=30); % indices of spikes
firings=[firings; t+0*fired,fired];
v(fired)=c(fired);
u(fired)=u(fired)+d(fired);
I=I+sum(S(:,fired),2);
v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
u=u+a.*(b.*v-u); % stability
end

%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
figure;
scatter(firings(:,1),firings(:,2),7,'filled'); 
xlabel('time (ms)');
ylabel('neuron');
title('Raster plot of activity');

%% SYNC

SIM_TIME=1000;
WINDOW=100;
num_active=zeros(1,10);
firings_temp=firings;
RASTER=zeros(Ne+Ni,SIM_TIME/WINDOW);
for j=1:10
    firings_temp=firings;
firings_temp(firings_temp(:,1)>WINDOW*j,:)=[];
firings_temp(firings_temp(:,1)<WINDOW*(j-1),:)=[];
neurons=firings_temp(:,2);
num_active(j)=length(unique(neurons))/(Ne+Ni);
end
sync=max(num_active);
%num_active=unique(firings_temp(:,2));
%for i=1:length(firings)
 %  RASTER(i,1)=firings(i,1); 
%end
%end

