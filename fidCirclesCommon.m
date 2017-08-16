function [ fids_common ] = fidCirclesCommon(density, radius, T2, omega, gammaG, t, TE)
%FIDCIRCLE Summary of this function goes here
%   density: the densities of circles
%   radius: the radii of circles in meter
%   T2: T2's of circles
%   omega: the frequency offsets of circles
%   gammaG = (gyromagnetic ratio) * (gradient strength)
%   t: time, a vector
% 

x = gammaG * t(:) * (radius(:)');

denrad = pi * density .* radius.^2;

integral = repmat( denrad(:)', [length(t), 1]) .* ( besselj(0, x) + besselj(2, x) );

phaseOmega = exp( -1i * t(:) * (omega(:)') );

decayT2 = exp( - ( t(:) + TE ) ./ (T2(:)') );

fids_common = integral .* phaseOmega .* decayT2;

end

