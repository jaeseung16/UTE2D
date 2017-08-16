function [ fid ] = sumOverCircles( location, density, radius, T2, omega, gammaG, t, nProj )
%SUMOVERCIRCLES Summary of this function goes here
%   location: the locations of circles in (r, theta)
%   density: the densities of circles
%   radius: the radii of circles in meter
%   T2: T2's of circles
%   omega: the frequency offsets of circles
%   gammaG = (gyromagnetic ratio) * (gradient strength)
%   t: time, a vector
% 

fids_common = fidCirclesCommon(density, radius, T2, omega, gammaG, t);

r0 = location(:,1);
theta0 = location(:,2);


kr0 = gammaG * t(:) * (r0(:)');

theta = theta0 + repmat( 2.0 * pi * ( 0:(nProj-1) ) / nProj, [length(theta0), 1] ) ;

for k = 1:length(r0)
    kr0_each = kr0(:,k);
    theta_each = theta(k,:);
    phase_each(:,k,:) = kr0_each * cos(theta_each);
end

fid = repmat( fids_common, [1, 1, nProj]);

phase = exp(1i * phase_each );

fid = phase .* fids_common;

fid = squeeze(sum(fid, 2));

end


