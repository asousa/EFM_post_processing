% Envelope peak method mostly taken from matlab's envelope function, except
% this version accepts calling findpeaks without minpeakdistance and gets a
% huge!!!! CPU time save.

function [yupper,ylower] = envPeak(x,iPkpos,iPkneg)

% pre-allocate space for results
nx = size(x,1);
yupper = zeros(size(x),'like',x);
ylower = zeros(size(x),'like',x);

% handle default case where not enough input is given
if nx < 2
  yupper = x;
  ylower = x;
  return
end

% compute upper envelope
for chan=1:size(x,2)
  if numel(iPkpos)<2
    % include the first and last points
    iLocs = [1; iPkpos; nx];
  else
    iLocs = iPkpos;
  end

  % smoothly connect the maxima via a spline.
  yupper(:,chan) = interp1(iLocs,x(iLocs,chan),(1:nx)','spline');
end

% compute lower envelope
for chan=1:size(x,2)
  if numel(iPkneg)<2
    % include the first and last points
    iLocs = [1; iPkneg; nx];
  else
    iLocs = iPkneg;
  end
  
  % smoothly connect the minima via a spline.
  ylower(:,chan) = interp1(iLocs,x(iLocs,chan),(1:nx)','spline');
end

