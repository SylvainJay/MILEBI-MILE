% _______________________________________________________________________
%
% assess_milebi_loglikelihood.m (August 23, 2017)
% _______________________________________________________________________

function L = assess_milebi_cost_function(Delta,rmeas,wl,Rb,model_spectra,angles,Gammas,Gammab)
% _______________________________________________________________________
% 
% This function computes the MILEBI cost function to be minimized.
% Inputs:   - Delta             : model input parameters,
%           - rmeas             : measured subsurface remote-sensing
%                                   reflectance
%           - wl                : wavebands,
%           - Rb                : bottom reflectances,
%           - model_spectra     : optical constants,
%           - angles            : acquisition geometry,
%           - Gammas            : environmental noise covariance matrix,
%           - Gammab            : bottom covariance matrices,
% Output:   - L                 : cost function value.
% _______________________________________________________________________
% 
% Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. 
% & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: 
% considering environmental noise and bottom intra-class variability for modeling 
% and inversion of water reflectance, Remote Sensing of Environment, 200(C), 352-367.
% _______________________________________________________________________
% 
B=Delta(5);
[rmod,Kb] = simu_nonnoisy_r(Delta,wl,Rb,model_spectra,angles);
Gamma_milebi = diag(Kb)*(B^2*Gammab(:,:,1)+(1-B)^2*Gammab(:,:,2))*diag(Kb) + Gammas;
L = sum(log(eig(Gamma_milebi))) + (rmeas-rmod)'*(Gamma_milebi\(rmeas-rmod));

end