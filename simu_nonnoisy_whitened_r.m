% _______________________________________________________________________
%
% simu_nonnoisy_whitened_r.m (August 23, 2017)
% _______________________________________________________________________

function whitened_rmod = simu_nonnoisy_whitened_r(Delta,wl,sqrt_inv_Gammas,Rb,model_spectra,angles)
% _______________________________________________________________________
% 
% This function simulates the subsurface remote-sensing reflectance of the
% water according to the bio-optical model used by Jay et al. (2017) and 
% originally introduced by Lee et al. (1998). This reflectance is then
% whitened using the environmental noise covariance matrix.
% Inputs:   - Delta             : model input parameters,
%           - wl                : wavebands,
%           - sqrt_inv_Gammas   : square root of the inverse of the
%                                   environmental noise covariance matrix,
%           - Rb                : bottom reflectances,
%           - model_spectra     : optical constants,
%           - angles            : acquisition geometry,
% Output:  - whitened_rmod     : whitened simulated reflectance,
% _______________________________________________________________________
% 
% Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. 
% & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: 
% considering environmental noise and bottom intra-class variability for modeling 
% and inversion of water reflectance, Remote Sensing of Environment, 200(C), 352-367.
% _______________________________________________________________________
% 
rmod = simu_nonnoisy_r(Delta,wl,Rb,model_spectra,angles);
whitened_rmod = sqrt_inv_Gammas*rmod;