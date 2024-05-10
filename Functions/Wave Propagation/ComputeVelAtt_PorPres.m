function [vp_pp, attp_pp] = ComputeVelAtt_PorPres(w, Material)
% Compute velocity and attenuation for compressional waves based on the
% porosity-pressure equations in dCS theory
% Reference: Spanos (2002), The Theormophysics of Porous Media.
% ------------------------------------------------------------------------

%% Compute polynomial constants
[AdCS, BdCS, CdCS] = getWaveCoeffs_dCSPorPres(Material, w);
% zero vector
Z = zeros(length(w),1);

%% Compute polynomial roots
pol_p = [AdCS' Z BdCS' Z CdCS'];
roots_p = zeros(length(w),4);
% loop over frequency values
for i = 1:length(w)
    % roots
    roots_p(i,:) = roots(pol_p(i,:));
end
% phase velocities
vp_pp = w'./real(roots_p);
% attenuations
attp_pp = -2*imag(roots_p)./real(roots_p);
