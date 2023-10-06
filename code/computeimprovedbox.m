function [improvedbox]= computeimprovedbox(controllerisfeasibility,channel,box,tol)
improvedlowerlimit = box(:,1);
improvedupperlimit = box(:,2);
for dim = 1:length(improvedlowerlimit)
    t = improvedlowerlimit;
    t(dim) = improvedupperlimit(dim);
    % if the temporary improvedupperlimit of the box  is not feasible, do bisection search to find the tigher uppper limit
    %if((~isfeasible(channel,power,Po,PAeff,t))&&(isfeasible(channel,power,Po,PAeff,improvedlowerlimit)))
    [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
    if(problem ~= 0) % infeasible problem
        tmax = improvedupperlimit(dim);
        tmin = improvedlowerlimit(dim);
        while (1)
            t(dim) = (tmax+tmin)/2;
            %if(isfeasible(channel,power,Po,PAeff,t))
            [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
            if(problem ==0) % if feasible
                tmin = t(dim);
                if((tmax-tmin) < tol)
                    break;
                end
            else
                tmax = t(dim);

            end
        end
    end
    improvedupperlimit(dim) = t(dim);
end

% Compress further
for dim = 1:length(improvedlowerlimit)
    t = improvedupperlimit;
    t(dim) = improvedlowerlimit(dim);
    %if((isfeasible(channel,power,Po,PAeff,t))&&(~isfeasible(channel,power,Po,PAeff,improvedupperlimit)))
    %if((isfeasible(channel,power,Po,PAeff,t)))
    [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
    if(problem == 0)
        tmax = improvedupperlimit(dim);
        tmin = improvedlowerlimit(dim);
        while (1)
            t(dim) = (tmax+tmin)/2;
            %if(isfeasible(channel,power,Po,PAeff,t))
            [~,problem] = controllerisfeasibility{{channel,[t(1); sqrt(exp(t(2:end))-1)]}};
            if(problem==0)
                tmin = t(dim);
                if((tmax-tmin) < tol)
                    break;
                end
            else
                tmax = t(dim);
            end
        end
    end
    improvedlowerlimit(dim) = t(dim);
end
improvedbox = [improvedlowerlimit improvedupperlimit];

