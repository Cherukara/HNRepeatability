function MaskThresholded = GradientBasedThreshold(FieldMap, Mask, alpha)
% Oliver C. Kiersnowski
% 1/2/2023
%
% Modified MT Cherukara (2023-02-14) to change the way the standard
%       deviation is calculated, assuming that the distribution of gradient
%       magnitudes follows a half-normal distribution (it actually doesn't,
%       it is an exponential, but the half-normal standard deviation is
%       nicer to deal with).

% Calculate the gradients in all three directions using finite difference
[Gx,Gy,Gz] = imgradientxyz(FieldMap);

% Calculate magnitude of gradient
Gmag = sqrt(Gx.^2 + Gy.^2 + Gz.^2);

% Work out standard deviation, with a correction for the half-normal distribution
StdGmag = std(Gmag(:))./(1-2/pi);

% Work out mean, also with a correction for the half-normal distribution
% (the value of the corrected mean should be very close to 0, so in theory
% this bit is not necessary)
MeanGmag = mean(Gmag(:)) - sqrt(2/pi).*StdGmag;

% Find the threshold
Thresh = MeanGmag + alpha*StdGmag;

% Create the mask
MaskThresholded = Mask;
MaskThresholded(Gmag > Thresh) = 0;

end