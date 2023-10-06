function [beamformers,currentbestUBstored,currentbestobjectivestored] = EEmaxMISOdownlinkoptimalBB(channel,sinrthreshold,power,Po,PAeff,tol)
[nUsers,nTx] = size(channel);
beamformers = zeros(nTx,nUsers);


tmin = [1/(Po+power/PAeff);log(1+sinrthreshold)];

tmax = [1/Po;log(1+power*diag(channel*channel'))];

newinitalbox = computeimprovedbox(channel,power,Po,PAeff,[tmin tmax],1e-3);
gap = 0;

ms = 3;

currentbestobj = newinitalbox(1,1)*sum(newinitalbox(2:end),1);
currentbestUB = newinitalbox(1,2)*sum(newinitalbox(2:end),2);

sequenceoflowerbounds = currentbestobj;
sequenceofupperbounds = currentbestUB;



currentbestobjectivestored = currentbestobj;
currentbestUBstored = currentbestUB;

nBoxes = 1;
boxstored(:,:,1) = newinitalbox;

% subplot(3,1,1)
% hold on
% initbox =  [tmin(2:end) tmax(2:end)];
% 
% set(0, 'defaultlinelinewidth',ms)
% rectangle('Position',[initbox(1,1)+gap initbox(2,1)+gap initbox(1,2)-initbox(1,1)-gap  initbox(2,2)-initbox(2,1)-gap],'edgecolor','r')
% xlim([initbox(1,1)-0.1 initbox(1,2)+0.1])
% ylim([initbox(2,1)-0.1 initbox(2,2)+0.1])
% xlabel('t1')
% ylabel('t2')
% subplot(3,1,2)
% hold on
% initbox =  [tmin(1:2) tmax(1:2)];
% 
% rectangle('Position',[initbox(1,1)+gap initbox(2,1)+gap initbox(1,2)-initbox(1,1)-gap  initbox(2,2)-initbox(2,1)-gap],'edgecolor','r')
% xlim([initbox(1,1) initbox(1,2)])
% ylim([initbox(2,1)-0.1 initbox(2,2)+0.1])
% xlabel('t0')
% ylabel('t1')
% 
% subplot(3,1,3)
% hold on
% initbox =  [tmin([1 3]) tmax([1 3])];
% 
% rectangle('Position',[initbox(1,1)+gap initbox(2,1)+gap initbox(1,2)-initbox(1,1)-gap  initbox(2,2)-initbox(2,1)-gap],'edgecolor','r')
% xlim([initbox(1,1) initbox(1,2)])
% ylim([initbox(2,1)-0.1 initbox(2,2)+0.1])
% xlabel('t0')
% ylabel('t2')

nColors = 20;
cc = hsv(nColors);
iColor = 1;
nIterations = 1;
while((currentbestUB - currentbestobj)> tol)
    currentbestUB - currentbestobj
    for iBox=1:nBoxes
%         iBox
%         sequenceoflowerbounds
%         nBoxes
%         currentbestUBstored
%         sequenceofupperbounds
        if(sequenceofupperbounds(iBox) == currentbestUB)
            % Branching
            lowerlimit = boxstored(:,1,iBox);
            upperlimit = boxstored(:,2,iBox);
            
            [temp,dim] = max(upperlimit-lowerlimit); % find the largest dimension for branching
            
            newupperlimit = upperlimit;
            newupperlimit(dim) = (upperlimit(dim)+lowerlimit(dim))/2; % calculate the upper limit for the new box
            
            newlowerlimit = lowerlimit;
            newlowerlimit(dim)=(upperlimit(dim)+lowerlimit(dim))/2; % calculate the upper limit for the new box
            
            newbox1 = [lowerlimit newupperlimit ]; % form the 1st box
            newbox2 = [newlowerlimit upperlimit ]; % form the 2nd box
            %% plot boxes
%             subplot(3,1,1)
%             box1 = newbox1(2:3,:);
%             box2 = newbox2(2:3,:);
%             %myplotbox(box1,'-ro',ms,0)
%             
%             
%             h11 = rectangle('Position',[box1(1,1)+gap  box1(2,1)+gap  box1(1,2)-box1(1,1)  box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'FaceColor',cc(mod(iColor,nColors/2)+1,:));
%             
%             plot(box1(1,1),box1(2,1),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             
%             
%             h12 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'FaceColor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             
% 
%             subplot(3,1,2)
%             box1 = newbox1(1:2,:);
%             box2 = newbox2(1:2,:);
%             
%             h21 = rectangle('Position',[box1(1,1)+gap  box1(2,1)+gap  box1(1,2)-box1(1,1)  box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'FaceColor',cc(mod(iColor,nColors/2)+1,:));
%             
%             plot(box1(1,1),box1(2,1),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             
%             
%             h22 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'FaceColor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'-s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'-s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             
% %             rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w')
% %             plot(box1(1,1)+gap,box1(2,1)+gap,'-ro','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% %             plot(box1(1,2),box1(2,2),'-ro','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% %             
% %             rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor','w')
% %             plot(box2(1,1)+gap,box2(2,1)+gap,'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% %             plot(box2(1,2),box2(2,2),'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%              
%             
%             subplot(3,1,3)
%             box1 = newbox1([1 3],:);
%             box2 = newbox2([1 3],:);
%             
%             h31 = rectangle('Position',[box1(1,1)+gap  box1(2,1)+gap  box1(1,2)-box1(1,1)  box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'FaceColor',cc(mod(iColor,nColors/2)+1,:));
%             
%             plot(box1(1,1),box1(2,1),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',3);
%             
%             %set(h1,'Visible','off');
%             
%             h32 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'FaceColor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'-s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'-s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
            
            
            %% clear boxes for the next plot
            
%             subplot(3,1,1)
%             box1 = newbox1(2:end,:);
%             box2 = newbox2(2:end,:);
%             set(h11,'Visible','off');
%             set(h12,'Visible','off');
%             set(h21,'Visible','off');
%             set(h22,'Visible','off');
%             set(h31,'Visible','off');
%             set(h32,'Visible','off');
% 
%             %rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w')
%             plot(box1(1,1)+gap,box1(2,1)+gap,'o','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             
%             %rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor','w')
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% 
%             subplot(3,1,2)
%             box1 = newbox1(1:2,:);
%             box2 = newbox2(1:2,:);
%             %rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w')
%             plot(box1(1,1)+gap,box1(2,1)+gap,'o','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             
%             %rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor','w')
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% 
%             subplot(3,1,3)
%             box1 = newbox1([1 3],:);
%             box2 = newbox2([1 3],:);
%             %rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w')
%             plot(box1(1,1)+gap,box1(2,1)+gap,'-ro','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box1(1,2),box1(2,2),'-ro','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             
%             %rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor','w')
%             plot(box2(1,1)+gap,box2(2,1)+gap,'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%             plot(box2(1,2),box2(2,2),'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);

            %% Improve the boxes
            
            improvednewbox1 = computeimprovedbox(channel,power,Po,PAeff,newbox1,1e-3);
            improvednewbox2 = computeimprovedbox(channel,power,Po,PAeff,newbox2,1e-3);
            
%             isfeasible(channel,power,Po,PAeff,improvednewbox1(:,1))
%             isfeasible(channel,power,Po,PAeff,improvednewbox1(:,1)+[0;0.02;0])
%             isfeasible(channel,power,Po,PAeff,improvednewbox2(:,1))
%             isfeasible(channel,power,Po,PAeff,improvednewbox2(:,1)+[0;0.02;0])
            %% plot improved boxes
%             subplot(3,1,1)
%             box1 = improvednewbox1(2:end,:);
%             box2 = improvednewbox2(2:end,:);
%             
%             h11 = rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'facecolor',cc(mod(iColor,nColors/2)+1,:));
%             plot(box1(1,1)+gap,box1(2,1)+gap,'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             
%             h12 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'facecolor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             
%             subplot(3,1,2)
%             box1 = improvednewbox1(1:2,:);
%             box2 = improvednewbox2(1:2,:);
%             
%             h21 = rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'facecolor',cc(mod(iColor,nColors/2)+1,:));
%             plot(box1(1,1)+gap,box1(2,1)+gap,'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             
%             h22 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'facecolor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             
%             subplot(3,1,3)
%             box1 = improvednewbox1([1 3],:);
%             box2 = improvednewbox2([1 3],:);
%             
%             h31 = rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+1,:),'facecolor',cc(mod(iColor,nColors/2)+1,:));
%             plot(box1(1,1)+gap,box1(2,1)+gap,'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             plot(box1(1,2),box1(2,2),'o','MarkerEdgeColor',cc(mod(iColor,nColors/2)+1,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+1,:),'Markersize',ms);
%             
%             h32 = rectangle('Position',[box2(1,1)+gap box2(2,1)+gap box2(1,2)-box2(1,1) box2(2,2)-box2(2,1)],'edgecolor',cc(mod(iColor,nColors/2)+11,:),'facecolor',cc(mod(iColor,nColors/2)+11,:));
%             plot(box2(1,1)+gap,box2(2,1)+gap,'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
%             plot(box2(1,2),box2(2,2),'s','MarkerEdgeColor',cc(mod(iColor,nColors/2)+11,:),'MarkerFaceColor',cc(mod(iColor,nColors/2)+11,:),'Markersize',ms);
            
            %% store boxes
            boxstored(:,:,iBox) = improvednewbox1;
            boxstored(:,:,1+nBoxes) = improvednewbox2;
            iColor = iColor+1;
            %% Bounding
            sequenceofupperbounds(iBox) = improvednewbox1(1,2)*sum(improvednewbox1(2:end,2));
            sequenceofupperbounds(1+nBoxes) = improvednewbox2(1,2)*sum(improvednewbox2(2:end,2));
            
            sequenceofupperboundscr = computeupperbound(channel,power,Po,PAeff,improvednewbox1);
            
            sequenceoflowerbounds(iBox) = improvednewbox1(1,1)*sum(improvednewbox1(2:end,1));
            sequenceoflowerbounds(1+nBoxes) = improvednewbox2(1,1)*sum(improvednewbox2(2:end,1));
            
            % Update global lower and upper bounds
            currentbestobj = max(sequenceoflowerbounds);
            currentbestUB = max(sequenceofupperbounds);
            
            % Increase the number of boxes
            nBoxes = nBoxes +1;
            
            %% Pruning
            prunedboxes = find(currentbestobj > sequenceofupperbounds);
            % clear the pruned boxes
            if(~isempty(prunedboxes))
                length(prunedboxes)
%             for n = 1:length(prunedboxes)
%                 boxtemp = boxstored(:,:,prunedboxes(n));
%                 box1 = boxtemp(2:end,:);
%                 subplot(3,1,1)
%                 rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w','facecolor','w')
%                 plot(box1(1,1)+gap,box1(2,1)+gap,'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%                 plot(box1(1,2),box1(2,2),'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%                 
%                 box1 = boxtemp(1:2,:);
%                 subplot(3,1,2)
%                 rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w','facecolor','w')
%                 plot(box1(1,1)+gap,box1(2,1)+gap,'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%                 plot(box1(1,2),box1(2,2),'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%                 box1 = boxtemp([1 3],:);
%                 subplot(3,1,3)
%                 rectangle('Position',[box1(1,1)+gap box1(2,1)+gap box1(1,2)-box1(1,1) box1(2,2)-box1(2,1)],'edgecolor','w','facecolor','w')
%                 plot(box1(1,1)+gap,box1(2,1)+gap,'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
%                 plot(box1(1,2),box1(2,2),'-ws','MarkerEdgeColor','w','MarkerFaceColor','w','Markersize',ms);
% 
%             end
            end
            
            boxstored(:,:,prunedboxes) = []; % delete the 
            sequenceoflowerbounds(prunedboxes) = [];
            sequenceofupperbounds(prunedboxes) = [];
            nBoxes = nBoxes-length(prunedboxes);
            nIterations = nIterations + 1;
            currentbestobjectivestored(nIterations) = currentbestobj;
            currentbestUBstored(nIterations) = currentbestUB;
            break;
        end
        
    end
    %save(savefile);
    
end
iBox = find(sequenceoflowerbounds==currentbestobj);
t = boxstored(:,1,iBox);
cvx_begin   quiet
variable    beamformers(nTx,nUsers) complex
subject to
for iUser=1:nUsers
    imag(channel(iUser,:)*beamformers(:,iUser)) == 0;
    real(channel(iUser,:)*beamformers(:,iUser)) >= ...
        sqrt(exp(t(iUser+1))-1)*norm([1 channel(iUser,:)*beamformers(:,1:nUsers~=iUser)]) ;
end
norm(vec(beamformers)) <= sqrt(min(power,PAeff*(1/t(1)-Po)));
cvx_end
cvx_status


