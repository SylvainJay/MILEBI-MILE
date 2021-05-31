% _______________________________________________________________________
% main.m
% version 1 (August 23, 2017)
% subroutines required: simu_nonnoisy_r.m, generate_lut_inputs.m
%                       latin_hs.m, ltqnorm.m, lhsu.m,
%                       simu_nonnoisy_whitened_r.m, assess_milebi_cost_function.m
% Matlab toolboxes required: Optimization Toolbox, Statistics Toolbox
% author: Sylvain Jay (sylvain.jay@fresnel.fr) 
% _______________________________________________________________________
% 
% This script runs MILE (MaxImum Likelihood estimation including
% Environmental noise) and MILEBI (MaxImum Likelihood estimation including 
% Environmental noise and Bottom Intra-class variability) (Jay et al., 2017) 
% in forward and inverse modes. 
% MILE and MILEBI enable one 
%       (1) to simulate realistic (i.e., as measured from a remote-sensing 
%       multi- or hyperspectral sensor) data sets of shallow water subsurface 
%       remote-sensing reflectance based on a given bio-optical model (e.g. 
%       the model of Lee et al. (1998)),
%       (2) to retrieve water column geophysical variables based on model 
%       inversion. In this script, these variables are the variables of Lee
%       et al.'s model (using a sum-to-one constrained linear mixture of two 
%       substrates to model the bottom reflectance), i.e.
%           - H     : Depth,
%           - P     : Phytoplankton absorption coefficient at 440 nm,
%           - G     : CDOM absorption coefficient at 440 nm,
%           - X     : Particle backscattering coefficient at 550 nm,
%           - B     : Bottom mixture coefficient.
% Note that MILE and MILEBI require the environmental noise matrix to be estimated
% beforehand. In addition, MILEBI requires the covariance matrix of
% each bottom class to be pre-estimated as well. Yet, for demonstration
% purposes, the matrices (and bottom mean reflectances) used by Jay et al. 
% (2017) are provided.
% For further details, please refer to the corresponding publication :
% 
% Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. 
% & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: 
% considering environmental noise and bottom intra-class variability for modeling 
% and inversion of water reflectance, Remote Sensing of Environment, 200(C):352-367.
%
% This work was supported by the French Defense Procurement Agency (DGA) with 
% the reference ANR-15-ASTR-0019 (HypFoM).
% _______________________________________________________________________

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load optical constants, pre-estimated parameters & define acquisition geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('wl.mat'); % Wavebands

load('water_optical_constants.mat'); % Water optical constants
% The data in water_optical_constants.mat are defined according to Jay et
% al. (2017) and the references provided therein:
%     -  model_spectra = [aw a0 a1 ag* bbw bbp*]
%           aw      : Pure water absorption coefficient (m-1) (Buiteveld 
%                       et al., 1994).
%           a0, a1  : Empirical spectra used to model phytoplankton 
%                       absorption (unitless) (Lee et al., 1998).
%           ag*     : CDOM specific absorption coefficient (unitless) 
%                       (e.g., Dekker et al., 2011).
%           bbw     : Pure water backscattering coefficient (m-1) (Morel, 1974).
%           bbp*    : Particle specific backscattering coefficient 
%                       (unitless) (Dekker et al., 2011; Lee et al., 2001).

load('Gammas.mat');  % Environmental noise covariance matrix
load('Rb.mat'); % Bottom mean reflectances (Sand, Oyster bags, Brown algae)
load('Gammab.mat'); % Bottom covariance matrices (Sand, Oyster bags, Brown algae)
figure,hold on,plot(wl,Rb(:,1),'LineWidth',2.5,'Color',[0.75,0.75,0.99]);plot(wl,Rb(:,2),'LineWidth',2.5,'Color',[1,.7,0]);plot(wl,Rb(:,3),'LineWidth',2.5,'Color',[0,1,.6]);legend('Sand','Oyster bags','Brown algae','Location','NorthWest');
plot(wl,mvnrnd(Rb(:,1)',squeeze(Gammab(:,:,1)),10),':','Color',[0,0,0.7]);plot(wl,mvnrnd(Rb(:,2)',squeeze(Gammab(:,:,2)),10),':','Color',[.75,.5,0]);plot(wl,mvnrnd(Rb(:,3)',squeeze(Gammab(:,:,3)),10),':','Color',[0,.75,.4]);
xlabel('Wavelength (nm)','FontWeight','Bold');ylabel('Remote-sensing reflectance','FontWeight','Bold'); title('Mean bottom reflectances & Simulated bottom intra-class variabilities','FontWeight','Bold'); 

angles = [30 0]; % [subsurface solar zenith angle / subsurface viewing zenith angle]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Simulations using MILE and MILEBI probabilistic forward models 
%           (sum-to-one constraint on bottom coefficients)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta = [2 0.1 0.1 0.01 .5]; % [H P G X B] : Define the values of bio-optical model parameters
Id_bott_1=1; Id_bott_2=3; % Define the indices of the considered bottoms (in Rb and Gammab)
[rmod,Kb] = simu_nonnoisy_r(Delta,wl,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles); % Simulate reflectance

Size_simu_image=4; Nb_simu=Size_simu_image*Size_simu_image; % Number of simulations
rmeas_mile = mvnrnd(rmod',Gammas,Nb_simu)'; % MILE probabilistic modeling
rmeas_milebi = mvnrnd(rmod', diag(Kb)*(Delta(5)^2*Gammab(:,:,Id_bott_1)...
    +(1-Delta(5))^2*Gammab(:,:,Id_bott_2))*diag(Kb)+Gammas,Nb_simu)';  % MILEBI probabilistic modeling

figure,subplot(1,2,1); plot(wl,rmod,'LineWidth',2.5,'Color',[0.75,0.75,0.99]); hold on; plot(wl,rmeas_mile,'Color',[0,0,0.7],'LineStyle',':'); hold off;
xlabel('Wavelength (nm)','FontWeight','Bold'); ylabel('Subsurface remote-sensing reflectance','FontWeight','Bold'); title('MILE probabilistic forward modeling','FontWeight','Bold'); legend('Noise-free reflectance','Noisy reflectance','Location','SouthWest');
subplot(1,2,2); plot(wl,rmod,'LineWidth',2.5,'Color',[0.75,0.75,0.99]); hold on; plot(wl,rmeas_milebi,'Color',[0,0,0.7],'LineStyle',':'); hold off;
xlabel('Wavelength (nm)','FontWeight','Bold'); ylabel('Subsurface remote-sensing reflectance','FontWeight','Bold'); title('MILEBI probabilistic forward modeling','FontWeight','Bold'); legend('Noise-free reflectance','Noisy reflectance','Location','SouthWest'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   MILE and MILEBI inversions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable setting for model inversion
% LB/UB = [ H   P       G       X       B]
lb = [      0   0       0       0       0]; % Lower bounds for optimization
ub = [      15  0.5     0.5     0.08    1]; % Upper bounds for optimization

options_mile = optimset('Algorithm','trust-region-reflective','TolFun',1e-10); % Optimization algorithm for MILE
options_milebi = optimset('Algorithm','interior-point','TolFun',1e-11); % Optimization algorithm for MILEBI
Nb_bott_types = size(Gammab,3); % Number of bottom types
Nb_bott_pairs = nchoosek(Nb_bott_types,2); % Number of possible bottom pairs
Nb_var = (length(lb)-1)+Nb_bott_types; % Number of variables to be estimated
sqrt_inv_Gammas = sqrtm(inv(Gammas)); % For MILE inversion
sol_mile=zeros(Nb_var,Nb_simu); sol_milebi=zeros(Nb_var,Nb_simu); % Initialize the vectors of estimated variables
N_mean_lut = 100; % Number of LUT entries to be averaged for initialization
percent_average_best_sol = 1+eps; % Determine the number of best bottom pairs to 
%                                   be averaged (= eps if only the best
%                                   pair is retained)

matlabpool(4); % If possible (and necessary), enables the optimization process 
               %   to be performed using several MATLAB workers

%% Generation of Look-Up Tables used for initialization
size_lut=10000;
P0_init = generate_lut_inputs(size_lut,lb,ub); % Generate LUT entries based on a Latin hypercube sampling scheme
parfor i=1:size_lut % Loop over all LUT entries 
    rinit_12(:,i) = simu_nonnoisy_r(P0_init(i,:),wl,Rb(:,[1 2]),model_spectra,angles); % Generate LUT entry for bottom pair #1
    rinit_13(:,i) = simu_nonnoisy_r(P0_init(i,:),wl,Rb(:,[1 3]),model_spectra,angles); % Generate LUT entry for bottom pair #2
    rinit_23(:,i) = simu_nonnoisy_r(P0_init(i,:),wl,Rb(:,[2 3]),model_spectra,angles); % Generate LUT entry for bottom pair #3
end

%% Model inversion using MILE and MILEBI
parfor i = 1:Nb_simu % Loop over all image pixels / reflectance spectra
%     rmeas = rmeas_mile(:,i); % input = subsurface remote-sensing reflectance simulated using MILE model
    rmeas = rmeas_milebi(:,i); % input = subsurface remote-sensing reflectance simulated using MILEBI model
    
    % Inversion for bottom pair #1
    Id_bott_1=1; Id_bott_2=2; Gammab_test = cat(3,Gammab(:,:,Id_bott_1),Gammab(:,:,Id_bott_2));
    MSE=mean((repmat(rmeas,[1 size_lut])-rinit_12).^2,1); [M,I]=sort(MSE); P0=mean(P0_init(I(1:N_mean_lut),:)); % Initialize using the "N_mean_lut" closest spectra in the LUT
    [Deltahat_mile_1,Resid_mile_1] = lsqcurvefit(@simu_nonnoisy_whitened_r,P0,wl,sqrt_inv_Gammas*rmeas,lb,ub,options_mile,sqrt_inv_Gammas,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles); % MILE inversion
    [Deltahat_milebi_1,Resid_milebi_1] = fmincon(@assess_milebi_cost_function,P0,[],[],[],[],lb,ub,[],options_milebi,rmeas,wl,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles,Gammas,Gammab_test); % MILEBI inversion
    
    % Inversion for bottom pair #1
    Id_bott_1=1; Id_bott_2=3; Gammab_test = cat(3,Gammab(:,:,Id_bott_1),Gammab(:,:,Id_bott_2));   
    MSE=mean((repmat(rmeas,[1 size_lut])-rinit_13).^2,1); [M,I]=sort(MSE); P0=mean(P0_init(I(1:N_mean_lut),:)); % Initialize using the "N_mean_lut" closest spectra in the LUT
    [Deltahat_mile_2,Resid_mile_2] = lsqcurvefit(@simu_nonnoisy_whitened_r,P0,wl,sqrt_inv_Gammas*rmeas,lb,ub,options_mile,sqrt_inv_Gammas,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles); % MILE inversion
    [Deltahat_milebi_2,Resid_milebi_2] = fmincon(@assess_milebi_cost_function,P0,[],[],[],[],lb,ub,[],options_milebi,rmeas,wl,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles,Gammas,Gammab_test); % MILEBI inversion

    % Inversion for bottom pair #1
    Id_bott_1=2; Id_bott_2=3; Gammab_test = cat(3,Gammab(:,:,Id_bott_1),Gammab(:,:,Id_bott_2));   
    MSE=mean((repmat(rmeas,[1 size_lut])-rinit_23).^2,1); [M,I]=sort(MSE); P0=mean(P0_init(I(1:N_mean_lut),:)); % Initialize using the "N_mean_lut" closest spectra in the LUT
    [Deltahat_mile_3,Resid_mile_3] = lsqcurvefit(@simu_nonnoisy_whitened_r,P0,wl,sqrt_inv_Gammas*rmeas,lb,ub,options_mile,sqrt_inv_Gammas,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles); % MILE inversion
    [Deltahat_milebi_3,Resid_milebi_3] = fmincon(@assess_milebi_cost_function,P0,[],[],[],[],lb,ub,[],options_milebi,rmeas,wl,Rb(:,[Id_bott_1 Id_bott_2]),model_spectra,angles,Gammas,Gammab_test); % MILEBI inversion
        
    % Group all vectors of estimated variables into matrices.
    Deltahat_mile = [[Deltahat_mile_1 (1-Deltahat_mile_1(5)) 0];[Deltahat_mile_2 0 (1-Deltahat_mile_2(5))];[Deltahat_mile_3(1:4) 0 Deltahat_mile_3(5) (1-Deltahat_mile_3(5))]];
    Deltahat_milebi = [[Deltahat_milebi_1 (1-Deltahat_milebi_1(5)) 0];[Deltahat_milebi_2 0 (1-Deltahat_milebi_2(5))];[Deltahat_milebi_3(1:4) 0 Deltahat_milebi_3(5) (1-Deltahat_milebi_3(5))]];

    % Sort the latter matrices according to the cost function values at the optima.
    Resid_mile = [Resid_mile_1;Resid_mile_2;Resid_mile_3];  Resid_milebi = [Resid_milebi_1;Resid_milebi_2;Resid_milebi_3];
    [Resid_mile_sorted,I_mile]=sort(Resid_mile);  [Resid_milebi_sorted,I_milebi]=sort(Resid_milebi);  
    Deltahat_mile=Deltahat_mile(I_mile,:);    Deltahat_milebi=Deltahat_milebi(I_milebi,:);

    % Average the vectors of estimated variables according to an increasing number of best bottom pairs
    averaged_Deltahat_mile = [mean(Deltahat_mile(1:1,:),1);mean(Deltahat_mile(1:2,:));mean(Deltahat_mile(1:3,:))];    
    averaged_Deltahat_milebi = [mean(Deltahat_milebi(1:1,:),1);mean(Deltahat_milebi(1:2,:));mean(Deltahat_milebi(1:3,:))];    
    
    % Determine the number of best bottom pairs to be averaged 
    Id_best_sol_mile = (abs((Resid_mile_sorted(:)-Resid_mile_sorted(1))/Resid_mile_sorted(1))*100) < percent_average_best_sol;
    Id_best_sol_milebi = (abs((Resid_milebi_sorted(:)-Resid_milebi_sorted(1))/Resid_milebi_sorted(1))*100) < percent_average_best_sol;
    
    % Select the optimum solutions according to the above-determined number of bottom pairs
    sol_mile(:,i)=averaged_Deltahat_mile(sum(Id_best_sol_mile),:)';
    sol_milebi(:,i)=averaged_Deltahat_milebi(sum(Id_best_sol_milebi),:)';    
end
matlabpool close;

