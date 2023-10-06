function energyefficiency = computeEE(channel,beamformer,Po,PAeff)
T = channel*beamformer;
signal = diag(T); % diagonal entries are h_k*w_k
signal = abs(signal).^2; % |h_k*w_k|^2 for each k
T(logical(eye(size(T)))) = 1; % replace the diagonals by 1
interference = sum(T.*conj(T),2); % norm of each row is the interference for each user
sumrate = sum(log(1 + signal./interference)); % rats in nat/s/Hz
energyefficiency = real(sumrate/(Po + (1/PAeff)*norm(beamformer(:))^2));