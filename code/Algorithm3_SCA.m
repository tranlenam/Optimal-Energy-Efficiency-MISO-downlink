function [EEbeamforming,beamformeropt,achievedpower,userrate]=Algorithm3_SCA(channel,sinrthreshold,power,Po,PAeff,beamformerinit)
[nUsers,nTx] = size(channel);
mybeta0 = zeros(nUsers,1);

myvartheta0 = zeros(nUsers,1);
myalphanext = zeros(nUsers,1);

for iUser = 1:nUsers
    mybeta0(iUser) = norm([1 channel(iUser,:)*beamformerinit(:,1:nUsers~=iUser)])^2;
    myvartheta0(iUser) = (real(channel(iUser,:)*beamformerinit(:,iUser))/sqrt(mybeta0(iUser)))^2+1 ;

    myalphanext(iUser) = log(myvartheta0(iUser));


end
z0 = (Po + (1/PAeff)*norm((beamformerinit(:)))^2)^2;
t0= sum(myalphanext)^2/z0;
nIterations = 20;
seqobj = zeros(nIterations+1,1);
EEbeamforming = zeros(nIterations+1,1);
EEbeamforming(1) = computeEE(channel,beamformerinit,Po,PAeff);

seqobj(1) = sqrt(t0);

beamformer = sdpvar(nTx,nUsers,'full','complex');
myvartheta = sdpvar(nUsers,1);
mybeta = sdpvar(nUsers,1); 
myalpha = sdpvar(nUsers,1);
sdpvar t z1 z
obj = t;
opts = sdpsettings('solver','mosek','verbose',0);
for iIteration =1:nIterations
    %
    F = [];
    for iUser = 1:nUsers
        % (8)
        F= [F,imag(channel(iUser,:)*beamformer(:,iUser)) == 0];

        %F = [F,real(channel(iUser,:)*beamformer(:,iUser)) >= sqrt(sinrthreshold(iUser))*...
        %             norm([1 channel(iUser,:)*beamformer(:,1:nUsers~=iUser)])];

        F = [F,cone([1 channel(iUser,:)*beamformer(:,1:nUsers~=iUser)],...
            real(channel(iUser,:)*beamformer(:,iUser))/sqrt(sinrthreshold(iUser)))];

        % (26b)
        F = [F,myvartheta(iUser) >= exp(myalpha(iUser))];

        %(27b)
        %F=[F,(mybeta(iUser)+1)/2 >= norm([1 channel(iUser,:)*beamformer(:,1:nUsers~=iUser) (mybeta(iUser)-1)/2])];
        F=[F,cone([1 channel(iUser,:)*beamformer(:,1:nUsers~=iUser) (mybeta(iUser)-1)/2],(mybeta(iUser)+1)/2 )];
        
        % (30b)
        %F=[F,real(channel(iUser,:)*beamformer(:,iUser)) >= (1/(2*mymu(iUser))*(myvartheta(iUser)-1)+(mymu(iUser)/2)*mybeta(iUser))];
        F=[F,real(channel(iUser,:)*beamformer(:,iUser)) >= firstorderapp(myvartheta(iUser)-1,mybeta(iUser),myvartheta0(iUser)-1,mybeta0(iUser))]; 
    end
    
    % (22c)
   
    F = [F,cone([z1;(z-1)/2],(z+1)/2)];
    
    F =[F,cone([sqrt(1/PAeff)*(beamformer(:));(z1-Po-1)/2],(z1-Po+1)/2)];
    
    F = [F,sum(myalpha) >= firstorderapp(t,z,t0,z0)];
    
    % (5b)
    F = [F,cone(beamformer(:),sqrt(power))];
    diagontics = optimize(F,-obj,opts); % solve the problem
    %
    if(diagontics.problem==0)||(diagontics.problem==4)
        beamformer0 = value(beamformer);
        t0 = value(t);
        z0 = value(z);
        myvartheta0 = value(myvartheta);
        mybeta0 = value(mybeta);
        seqobj(iIteration+1) = sqrt(value(obj));
        EEbeamforming (iIteration+1) = computeEE(channel,beamformer0,Po,PAeff);
        achievedpower = norm((beamformer0(:)))^2;
        userrate = zeros(nUsers,1);
        for iUser = 1:nUsers
            userrate(iUser) = (log(1+abs(channel(iUser,:)*beamformer0(:,iUser))^2 ...
                /(norm([1 channel(iUser,:)*beamformer0(:,1:nUsers~=iUser)])^2)));

        end
    else
        EEbeamforming = NaN; achievedpower = NaN; seqobj = NaN; beamformer0 = NaN;
    end
    %}

end
beamformeropt = beamformer0;