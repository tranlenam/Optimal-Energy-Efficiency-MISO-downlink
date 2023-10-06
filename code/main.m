clear 
clc
rng('default')
nUsers = 3; % number of users
nTx = 4; % number of tx antennas

%% power consumption
power = 46; % dB
power = 10.^((power-30)/10); % to linear scale
PAeff = 0.38; %  power amplifier efficiency
Pdyn = 42; % dynamic power consumption (dBm)
Pdyn = 10^((Pdyn-30)/10); % to linear scale
Psta = 33; % dBm
Psta = 10^((Psta-30)/10); % to linear scale
Po = nTx*Pdyn+Psta; % total circuit power


sinrthreshold = ones(nUsers,1); % SINR requirements (i.e., \bar{\gamma}_k in (4b)

%% path loss
d = 0.1+0.7*rand(nUsers,1) ; % users's distance: .1km to .8 km
PL = 128.1+37.6*log10(d); % pass loss model (in dB)
PL = 10.^(-PL/10);
B = 10e6; % bandwidth
No = -203.975; % dBW/Hz
noise_p = No+10*log10(B)+30; % (~104 dBm)

noise_p = 10^((noise_p-30)/10); % to linear scale

%% main loop here

% generate small scale fading channel
channel = sqrt(PL).*(sqrt(1/2)*(randn(nUsers,nTx)+1i*randn(nUsers,nTx)));
channel = channel/sqrt(noise_p); % normalize the channel with noise power

% the noise power in SINR expressions now becomes 1

%% SPmin power
[minpowerSPmin,beamformerSPmin]= SPmin(channel,sinrthreshold);
EEofSPmin = computeEE(channel,beamformerSPmin,Po,PAeff);

%% Algorithm 3 (SCA)

[EE_SCA_seq,beamformerEEmax_SCA]= Algorithm3_SCA(channel,sinrthreshold,power,Po,PAeff,beamformerSPmin);

%% Algorithm 1 (BRB)
% In Algorithm 1, we need to REPEATEDLY check if (10) is feasible or not.
% For this purpose, we first create an optimizer for (10)

channelpara = sdpvar(nUsers, nTx,'full','complex'); % treat channel in (10) as optimization variable 
ops = sdpsettings('solver','mosek','verbose',0);
beamformer = sdpvar(nTx,nUsers,'full','complex');
u = sdpvar(nUsers,nUsers,'full','complex');
t = sdpvar(nUsers+1,1); % also treat t as opt variable
myalpha = sdpvar(nUsers,1);

F = [u == channelpara*beamformer];

for iUser=1:nUsers
    F = [F,imag(u(iUser,iUser)) == 0]; % this is (10c)
%     F = [F,real(u(iUser,iUser)) >= t(iUser+1)*myalpha(iUser)]; % (10b) 
%     F = [F,cone([1 u(iUser,1:nUsers~=iUser)],myalpha(iUser))]; % (10b)
    F = [F,cone([1 u(iUser,1:nUsers~=iUser)],real(u(iUser,iUser))/t(iUser+1))]; % (10b)
end
F = [F,cone(beamformer(:),sqrt(min(power,PAeff*(1/t(1)-Po))))]; % (10d)

controllerisfeasibility = optimizer(F,[],ops,{channelpara,t},beamformer); % indicate channelpara and t are parameters of the optimizer
%}

tol = 0.005; % error tolerance

[optimalbeamformerEEmax,globalUB,globalLB] = Algorithm1_BRB(controllerisfeasibility,channel,sinrthreshold,power,Po,PAeff,tol);
computeEE(channel,optimalbeamformerEEmax,Po,PAeff)


semilogx(EE_SCA_seq*B/1e6*log2(exp(1)),'b')
hold on
semilogx(globalLB*B/1e6*log2(exp(1)),'r')
semilogx(globalUB*B/1e6*log2(exp(1)),'b')
ylabel('Energy Efficiency(Mb/J)')
xlabel('Iteration count')
legend('Alg. 3 (SCA)','Alg. 1 (BRB - Lower bound)','Alg. 1 (BRB - Upper bound)')
saveas(gcf,'../results/convergence.png')