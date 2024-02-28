%%  run file for Monte Carlo solution to radiative transfer equations
clear variables;close all;clc;
rng(645651)
possiblePaths = {'/homes/male7736/Desktop/Research/Monte-Carlo-Code',...
    'C:\Users\Test\Documents\Research\Monte-Carlo-Code'} ;
addpath( possiblePaths{ispc+1} );%Monte Carlo Code
input=struct();
%% Input parameters
ssa = 1;%Single scattering albedo: P(scattered)
nPhotons = 10000; %Number of photons to simulate
opticalDepth = 10; %Optical depth of the cloud
meanFreePath = 1; %Mean free path: Mean distance traveled before interaction 
g = 0;%Asymmetry parameter = E[cos(\theta)] where theta is the scattering angle
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

if strcmpi(model, 'twostream')
results = TwoStreamSimulate(input);
else
    error('Only the two stream model is implemented') 
end



%No absorption 
% g = input.scatteringProbs(1)-input.scatteringProbs(2);
tauBar = input.opticalDepth/input.mfp ;
Rtheory= tauBar*(1-g)*0.5/(1+tauBar*(1-g)*0.5);
Ttheory= 1/(1+tauBar*(1-g)*0.5);

t = 0:0.01:opticalDepth;
figure;plot(t, t*(1-g)*0.5./(1+t*(1-g)*0.5))
hold on;
plot(t, 1./(1+t*(1-g)*0.5))
legend('Reflectance',  'Transmittance')
xlabel('Optical depth')
ylabel('R or T')
figure;
plot( linspace(0,opticalDepth, size(results.FDwnUp,2)),results.FDwnUp )
xlabel('Depth');ylabel('F/F_0')
legend('F_{dwn}', 'F_{up}')






