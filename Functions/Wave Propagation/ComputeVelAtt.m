function [vp, attp, vs, atts] = ComputeVelAtt(w, Material)
% Compute velocity and attenuation for compressional and shear waves using
% different porous media models
% ------------------------------------------------------------------------

%% Compute polynomial constants
% Biot (BT) theory
[Ap_BT, Bp_BT, Cp_BT, As_BT, Bs_BT, Cs_BT] = getWaveCoeffs_BT(Material, w);
% Zhao - simplified de la Cruz and Spanos - (ZH) theory
[Ap_ZH, Bp_ZH, Cp_ZH, As_ZH, Bs_ZH, Cs_ZH] = getWaveCoeffs_ZH(Material, w);
% de la Cruz and Spanos (dCS) theory
[Ap_dCS, Bp_dCS, Cp_dCS, As_dCS, Bs_dCS, Cs_dCS] = getWaveCoeffs_dCS(Material, w);
% zero vector
Z = zeros(length(w),1);

% regroup constants - P wave
Ap = [Ap_BT; Ap_ZH; Ap_dCS];
Bp = [Bp_BT; Bp_ZH; Bp_dCS];
Cp = [Cp_BT; Cp_ZH; Cp_dCS];
% regroup constants - S wave
As = [As_BT; As_ZH; As_dCS];
Bs = [Bs_BT; Bs_ZH; Bs_dCS];
Cs = [Cs_BT; Cs_ZH; Cs_dCS];

%% Initialize variables
vp = zeros(length(w),6);
vs = zeros(length(w),6);
attp = zeros(length(w),6);
atts = zeros(length(w),6);

%% Compute polynomial roots
% loop over models
for m = 1:3
    % polynomial 4th order P wave
    pol_p = [Ap(m,:)' Z Bp(m,:)' Z Cp(m,:)'];
    % initialize matrices
    roots_p = zeros(length(w),4);
    kp = zeros(length(w),2);
    % loop over frequency values
    for i = 1:length(w)
        % roots
        roots_p(i,:) = roots(pol_p(i,:));
        % physical solutions
        kp(i,:) = roots_p(i, real(roots_p(i,:))>0);
    end
    % phase velocities
    vp(:, 2*m-1:2*m) = w'./real(kp);
    % attenuations
    attp(:, 2*m-1:2*m) = -2*imag(kp)./real(kp);
%     attp(:, 2*m-1:2*m) = -imag(kp);
    
    % polynomial 4th order S wave
    pol_s = [As(m,:)' Z Bs(m,:)' Z Cs(m,:)'];
    % initialize matrices
    roots_s = zeros(length(w),4);
    ks = zeros(length(w),2);
    % loop over frequency values
    for i = 1:length(w)
        % roots
        aux = roots(pol_s(i,:));
        roots_s(i,1:length(aux)) = aux;
        % physical solutions
        ks(i,1:length(aux)/2) = roots_s(i, real(roots_s(i,:))>0);
    end
    % phase velocities
    vs(:, 2*m-1:2*m) = w'./real(ks);
    % attenuations
    atts(:, 2*m-1:2*m) = -2*imag(ks)./real(ks);
%     atts(:, 2*m-1:2*m) = -imag(ks);
end

end