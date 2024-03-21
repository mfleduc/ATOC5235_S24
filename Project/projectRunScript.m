%%  run file for Monte Carlo solution to radiative transfer equations
clear variables;close all;clc;
rng(645651)
possiblePaths = {'/homes/male7736/Desktop/Research/Monte-Carlo-Code',...
    'C:\Users\Test\Documents\Research\Monte-Carlo-Code'} ;
addpath( possiblePaths{ispc+1} );%Monte Carlo Code
input=struct();
%% Input parameters
ssa = 1;%Single scattering albedo: P(scattered)
nPhotons = 8000; %Number of photons to simulate
opticalDepths = [0.001,0.5:0.5:20]; %Optical depths of the cloud
meanFreePath = 1; %Mean free path: Mean distance traveled before interaction 
g = 0;%Asymmetry parameter = E[cos(\theta)] where theta is the scattering angle
%Flags
model = 'twostream'; %twostream: 2 stream scattering model. multi: requires full scattering pattern
medIsIsotropic = true;
%% Compile the inputs
input.Nphotons = nPhotons;
input.ssa = ssa;

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

if strcmpi(model, 'twostream')
    for kk = 1:length(opticalDepths)
        input.opticalDepth = opticalDepths(kk);
        results = TwoStreamSimulate(input);
        rCalc(kk) = results.RTA(1);
        tCalc(kk) = results.RTA(2);
    end
else
    error('Only the two stream model is implemented') 
end
%No absorption 
Rtheory = opticalDepths*(1-g)*0.5./(1+opticalDepths*(1-g)*0.5) ;
Ttheory= 1./(1+opticalDepths*(1-g)*0.5);

figure;plot(opticalDepths,Rtheory )
hold on;
plot(opticalDepths,Ttheory )
plot(opticalDepths, rCalc,'b.')
plot(opticalDepths, tCalc, 'r.')
legend('Reflectance (theory)',  'Transmittance (theory)',...
    'Reflectance (Calculated)','Transmittance (Calculated)')
xlabel('Optical depth')
ylabel('R or T')
grid on
title(sprintf('Theoretical vs calculated reflectance/transmittance,g=%.2f',g))







