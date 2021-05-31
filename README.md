# MILEBI and MILE methods for hyperspectral remote sensing of shallow waters

These MATLAB scripts/functions run MILE (MaxImum Likelihood estimation including Environmental noise) and MILEBI (MaxImum Likelihood estimation including Environmental noise and Bottom Intra-class variability) (Jay et al., 2017) in forward and inverse modes.

MILE and MILEBI make it possible :
(1) to simulate realistic (i.e., as measured from a remote-sensing multi- or hyperspectral sensor) data sets of shallow water subsurface remote-sensing reflectance based on a given bio-optical model (e.g. the model of Lee et al. (1998)),
(2) to retrieve water column geophysical variables based on model inversion. In this script, these variables are the variables of Lee et al.'s model (using a sum-to-one constrained linear mixture of two substrates to model the bottom reflectance).

Note that MILE and MILEBI require the environmental noise matrix to be estimated beforehand. In addition, MILEBI requires the covariance matrix of each bottom class to be pre-estimated as well. Yet, for demonstration purposes, the matrices (and bottom mean reflectances) used by Jay et al. (2017) are provided.

For further details, please refer to the corresponding publication : 
Jay, S., Guillaume, M., Minghelli, A., Deville, Y., Chami, M., Lafrance, B. & Serfaty, V. (2017), Hyperspectral remote sensing of shallow waters: Considering environmental noise and bottom intra-class variability for modeling and inversion of water reflectance, Remote Sensing of Environment, 200(C):352-367.
