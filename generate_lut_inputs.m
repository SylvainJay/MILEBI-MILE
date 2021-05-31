% _______________________________________________________________________
%
% generate_lut_inputs.m (August 23, 2017)
% _______________________________________________________________________

function lut_inputs = generate_lut_inputs(size_lut,lb,ub)
% _______________________________________________________________________
% 
% This function generates numerous sets of LUT inputs parameters according
% the Latin Hypercube sampling scheme used by Jay et al. (2017) and originally 
% proposed by Garcia et al. (2014).
% Inputs:   - size_lut          : Number of LUT entries,
%           - lb                : Lower bounds,
%           - ub                : Upper bounds,
% Outputs:  - lut_inputs        : LUT inputs variables.
% _______________________________________________________________________
%
% Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. 
% & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: 
% considering environmental noise and bottom intra-class variability for modeling 
% and inversion of water reflectance, Remote Sensing of Environment, 200(C), 352-367.
%
% Garcia, R. A., McKinna, L. I., Hedley, J. D., & Fearns, P. R. (2014). 
% Improving the optimization solution for a semi?analytical shallow water 
% inversion model in the presence of spectrally correlated noise. Limnology 
% and Oceanography: Methods, 12(10), 651-669.
% _______________________________________________________________________
% 
% Variable setting
size_lut2=size_lut*20;
Nb_var = length(lb);

% Computation of means and standard deviations of Gaussian distributions 
% used to generate the H, P, G and X values.
step=1/500;
Y = normpdf((1/3),0,0:step:1)./normpdf(0,0,0:step:1);
step_h=(ub(1)-lb(1))*step; step_p=(ub(2)-lb(2))*step; step_g=(ub(3)-lb(3))*step; step_x=(ub(4)-lb(4))*step;
[m,Id]=min((Y-0.5).^2); 
std_h=step_h*(Id-1)+lb(1); std_p=step_p*(Id-1)+lb(2); std_g=step_g*(Id-1)+lb(3); std_x=step_x*(Id-1)+lb(4);
xmean = zeros(1,Nb_var-1); xsd = [std_h std_p std_g std_x];

% Generate LUT entries based on Latin Hypercube sampling (Gaussian
% distributions for H, P, G and X, and uniform distribution for B)
lut_inputs=latin_hs(xmean,xsd,size_lut2,4);
lut_inputs = lut_inputs( lut_inputs(:,1)>0 & lut_inputs(:,2)>0 & lut_inputs(:,3)>0 & lut_inputs(:,4)>0,:); 
Perm = randperm(size(lut_inputs,1),size_lut); 
lut_inputs=lut_inputs(Perm,:);
lut_inputs = [lut_inputs lhsu(lb(5),ub(5),size_lut)];

