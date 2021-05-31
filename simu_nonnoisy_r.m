% _______________________________________________________________________
%
% simu_nonnoisy_r.m (August 23, 2017)
% _______________________________________________________________________

function [rmod,Kb] = simu_nonnoisy_r(Delta,wl,Rb,model_spectra,angles)
% _______________________________________________________________________
% 
% This function simulates the subsurface remote-sensing reflectance of the
% water according to the bio-optical model used by Jay et al. (2017) and 
% originally introduced by Lee et al. (1998).
% Inputs:   - Delta             : model input parameters,
%           - wl                : wavebands,
%           - Rb                : bottom reflectances,
%           - model_spectra     : optical constants,
%           - angles            : acquisition geometry,
% Outputs:  - rmod              : simulated reflectance,
%           - Kb                : simulated Kb.
% _______________________________________________________________________
% 
% Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. 
% & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: 
% considering environmental noise and bottom intra-class variability for modeling 
% and inversion of water reflectance, Remote Sensing of Environment, 200(C), 352-367.
%
% Lee, Z., Carder, K. L., Mobley, C. D., Steward, R. G., & Patch, J. S. (1998). 
% Hyperspectral remote sensing for shallow waters. I. A semianalytical model. 
% Applied optics, 37(27), 6329-6338.
% _______________________________________________________________________
% 

H=Delta(1); P=Delta(2); G=Delta(3); X=Delta(4); B=Delta(5); Nb_wl=length(wl);

% Modeling of absorption and backscattering coefficients
aw = model_spectra(:,1);
aphy = (model_spectra(:,2)+model_spectra(:,3).*log(repmat(P,[Nb_wl 1]))).*repmat(P,[Nb_wl 1]);
acdom = repmat(G,[Nb_wl 1]) .* model_spectra(:,4);
a = aw + aphy + acdom;
bbw = model_spectra(:,5);
bbp = repmat(X,[Nb_wl 1]) .* model_spectra(:,6);
bb = bbw + bbp;

% Modeling of subsurface remote-sensing reflectance
k = (a+bb); 
x = bb./k;

rinf = (0.084+0.17*x) .* x;
kd = k/cosd(angles(1));
kub = (1/cosd(angles(2))) * (1.04*k) .* (1+5.4*x).^0.5;
kuc = (1/cosd(angles(2))) * (1.03*k) .* (1+2.4*x).^0.5;

Kb = exp(-(kd+kub).*repmat((H)',[Nb_wl 1]));
Kc = exp(-(kd+kuc).*repmat((H)',[Nb_wl 1]));

rmod = (ones(Nb_wl,1)-Kc).*rinf + Kb.*(repmat(B',[Nb_wl 1]).*Rb(:,1)+repmat((1-B)',[Nb_wl 1]).*Rb(:,2));
