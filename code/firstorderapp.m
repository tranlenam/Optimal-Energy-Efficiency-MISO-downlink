function output = firstorderapp(x,y,x0,y0)
% approximate \sqrt(xy) around (x0,y0)
%   Detailed explanation goes here
output = sqrt(x0*y0)+1/2*sqrt(x0/y0)*(y-y0)+1/2*sqrt(y0/x0)*(x-x0);
end