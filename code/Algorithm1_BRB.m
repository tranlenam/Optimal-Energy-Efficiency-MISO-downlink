function [optbeamformers,currentbestUBstored,currentbestobjectivestored] = Algorithm1_BRB(controllerisfeasibility,channel,sinrthreshold,power,Po,PAeff,tol)
[nUsers,nTx] = size(channel);

tmin = [1/(Po+power/PAeff);log(1+sinrthreshold)];

tmax = [1/Po;log(1+power*diag(channel*channel'))];

newinitalbox = computeimprovedbox(controllerisfeasibility,channel,[tmin tmax],1e-3);

currentbestobj = newinitalbox(1,1)*sum(newinitalbox(2:end,1));
currentbestUB = newinitalbox(1,2)*sum(newinitalbox(2:end,2));

sequenceoflowerbounds = currentbestobj;
sequenceofupperbounds = currentbestUB;
sequenceofimprovedupperbounds = currentbestUB;



currentbestobjectivestored = currentbestobj;
currentbestUBstored = currentbestUB;

nBoxes = 1;
boxstored(:,:,1) = newinitalbox;

maxIteration = 1e4;
for iIter =1:maxIteration
    [~,iBoxes] = max(sequenceofupperbounds);
    iBox = iBoxes(randi(length(iBoxes)));
    % Branching
    lowerlimit = boxstored(:,1,iBox);
    upperlimit = boxstored(:,2,iBox);

    [~,dim] = max(upperlimit-lowerlimit); % find the largest dimension for branching

    newupperlimit = upperlimit;
    newupperlimit(dim) = (upperlimit(dim)+lowerlimit(dim))/2; % calculate the upper limit for the new box

    newlowerlimit = lowerlimit;
    newlowerlimit(dim)=(upperlimit(dim)+lowerlimit(dim))/2; % calculate the upper limit for the new box

    newbox1 = [lowerlimit newupperlimit ]; % form the 1st box
    newbox2 = [newlowerlimit upperlimit ]; % form the 2nd box


    %% Improve the boxes
    t = newbox1(:,1);
    [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
    if(problem==0)
        improvednewbox1 = computeimprovedbox(controllerisfeasibility,channel,newbox1,1e-3);
    end

    t = newbox2(:,1);
    [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
    if(problem==0)
        improvednewbox2 = computeimprovedbox(controllerisfeasibility,channel,newbox2,1e-3);
    end


    %% store boxes
    boxstored(:,:,iBox) = improvednewbox1;
    boxstored(:,:,1+nBoxes) = improvednewbox2;
    % iColor = iColor+1;
    %% Bounding
    sequenceofupperbounds(iBox) = improvednewbox1(1,2)*sum(improvednewbox1(2:end,2));
    sequenceofupperbounds(1+nBoxes) = improvednewbox2(1,2)*sum(improvednewbox2(2:end,2));

    sequenceoflowerbounds(iBox) = improvednewbox1(1,1)*sum(improvednewbox1(2:end,1));
    sequenceoflowerbounds(1+nBoxes) = improvednewbox2(1,1)*sum(improvednewbox2(2:end,1));

    % Update global lower and upper bounds
    currentbestobj = max(sequenceoflowerbounds);
    currentbestUB = max(sequenceofupperbounds);

    % Increase the number of boxes
    nBoxes = nBoxes +1;

    %% Pruning
    prunedboxes = find(currentbestobj > sequenceofupperbounds);

    boxstored(:,:,prunedboxes) = []; % delete the
    sequenceoflowerbounds(prunedboxes) = [];
    sequenceofupperbounds(prunedboxes) = [];
    nBoxes = nBoxes-length(prunedboxes);
    currentbestobjectivestored(iIter+1) = currentbestobj;
    currentbestUBstored(iIter+1) = currentbestUB;
    if(abs(currentbestUB - currentbestobj) < tol)
        break
    end

end

iBox = find(sequenceoflowerbounds==currentbestobj);
t = boxstored(:,1,iBox);
opts  = sdpsettings('solver','mosek','verbose',0);
beamformers = sdpvar(nTx,nUsers,'full','complex');
F = [];
for iUser=1:nUsers
    F = [F,imag(channel(iUser,:)*beamformers(:,iUser)) == 0];
    F = [F,cone([1 channel(iUser,:)*beamformers(:,1:nUsers~=iUser)],...
        real(channel(iUser,:)*beamformers(:,iUser))/sqrt(exp(t(iUser+1))-1))];
end
F = [F,cone(beamformers(:),sqrt(min(power,PAeff*(1/t(1)-Po))))];
diagnotics = optimize(F,[],opts);
if (diagnotics.problem==0)||(diagnotics.problem==4)
    optbeamformers = value(beamformers);
else
    optbeamformers = NaN;
end
end


