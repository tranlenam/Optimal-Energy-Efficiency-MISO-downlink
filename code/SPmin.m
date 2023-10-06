function [minpower,beamformer]=SPmin(channel,sinrthreshold)

[nUsers,nTx] = size(channel);

ops = sdpsettings('solver','mosek','verbose',0); % use MOSEK
beamformer = sdpvar(nTx,nUsers,'full','complex'); % define optimization variable
F=[];
obj = norm(vec(beamformer));
for iUser=1:nUsers
    F = [F,imag(channel(iUser,:)*beamformer(:,iUser)) == 0;];
    F = [F,real(channel(iUser,:)*beamformer(:,iUser)) >=...
        sinrthreshold(iUser)...
        *norm([1 channel(iUser,:)*beamformer(:,1:nUsers~=iUser)])];
end

diagnotics = optimize(F,obj,ops);
if(diagnotics.problem==0)
    minpower = double(obj)^2;
    beamformer = value(beamformer);
else
    minpower = NaN; 
    beamformer = NaN;
end