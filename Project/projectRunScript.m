%%  run file for Monte Carlo solution to radiative transfer equations
clear variables;close all;clc;
rng(645651)
addpath( '/homes/male7736/Desktop/Research/Monte-Carlo-Code' );%Monte Carlo Code
input=struct();
%% Input parameters
ssa = 0.8;%Single scattering albedo: P(scattered)
nPhotons = 10000; %Number of photons to simulate
opticalDepth = 1; %Optical depth of the cloud
meanFreePath = 1; %Mean free path: Mean distance traveled before interaction 
g = -0.75;        %Asymmetry parameter = E[cos(\theta)] where theta is the scattering angle
%Flags
model = 'twostream'; %twostream: 2 stream scattering model. multi: requires full scattering pattern
medIsIsotropic = true;
%% Compile the inputs
input.Nphotons = nPhotons;
input.ssa = ssa;
input.opticalDepth = opticalDepth;
input.mfp = meanFreePath;
input.scatteringProbs = [1+g,1-g]/2; %If isotropic: [P(forward),P(backward)]
                                   %If not: [[P(forward|going down),P(backward|going down)];
                                   %         [P(forward|going up)  ,P(backward|going up)  ]]
                                   %In each row the 2nd entry is uniquely
                                   %determined by the first, and only the
                                   %first is used, but both should be
                                   %included for clarity.
%Flags
input.flags = struct;
input.flags.model = model;
input.flags.isotropicMedium = medIsIsotropic;
results = Simulate(input);

%No absorption
% g = input.scatteringProbs(1)-input.scatteringProbs(2);
Rtheory= opticalDepth*(1-g)*0.5/(1+opticalDepth*(1-g)*0.5);
Ttheory= 1/(1+opticalDepth*(1-g)*0.5);