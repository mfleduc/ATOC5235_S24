%%  run file for Monte Carlo solution to radiative transfer equations
clear variables;close all;clc;
rng(645651)
possiblePaths = {'/homes/male7736/Desktop/Research/Monte-Carlo-Code',...
    'C:\Users\Test\Documents\Research\Monte-Carlo-Code'} ;
addpath( possiblePaths{ispc+1} );%Monte Carlo Code
input=struct();
%% Input parameters
ssa = 1;%Single scattering albedo: P(scattered)
nPhotons = 100000; %Number of photons to simulate
opticalDepths = [ 1 ]; %Optical depths of the cloud
meanFreePath = 1; %Mean free path: Mean distance traveled before interaction 
g = 0;%Asymmetry parameter = E[cos(\theta)] where theta is the scattering angle
%Flags
model = 'twostream'; %twostream: 2 stream scattering model. 3d: requires full scattering pattern
medIsIsotropic = true;
%% Compile the inputs
input.Nphotons = nPhotons;
input.ssa = ssa;
input.opticalDepth = opticalDepths ;
input.mfp = meanFreePath;
input.scatteringProbs = 'Henyey-Greenstein';
input.asymmetry = g;
%Flags
input.flags = struct;
input.flags.model = model;
input.flags.isotropicMedium = medIsIsotropic;
if strcmpi(model, 'twostream')
    input.scatteringProbs = [1+g,1-g]/2; %If isotropic: [P(forward),P(backward)]
                                   %If not: [[P(forward|going down),P(backward|going down)];
                                   %         [P(forward|going up)  ,P(backward|going up)  ]]
                                   %In each row the 2nd entry is uniquely
                                   %determined by the first, and only the
                                   %first is used, but both should be
                                   %included for clarity.
    for kk = 1:length(opticalDepths)
        input.opticalDepth = opticalDepths(kk);
        results = TwoStreamSimulate(input);
        rCalc(kk) = nnz(results.nSteps==3&results.displacements<0)/nnz(results.nSteps>2) ;
        tCalc(kk) = nnz(results.nSteps==3&results.displacements>0)/nnz(results.nSteps>2) ;
    end
else
    input.scatteringProbs = 'Henley-Greenstein'; %Either the name of the 
    %phase function (only takes Henley-Greenstein) or a struct with fields
    %x,y,type={'pdf' or 'cdf'} that has x = scattering angles, y = pdf/cdf(x)
    for kk = 1:length(opticalDepths)
        input.opticalDepth = opticalDepths(kk);
        input.Nphotons = ceil(nPhotons*input.opticalDepth^3);
        results = ThreeDSimulate(input);
        rCalc(kk) = results.RTA(1);
        tCalc(kk) = results.RTA(2);
    end
end
%No absorption 
% %pr3steps = @(x,z)0.25*((2*x+1).*exp(-2*x)-exp(-2*z))./((1-exp(-z))*(2-exp(-x)-exp(x-z)));
% %pt3steps = @(x,z)0.25*exp(-z)/(1-exp(-z))*(2*(z-x)-exp(-2*x))./(2-exp(-x)-exp(x-z)) ;
% % pr3steps = @(x,z)(-1*(1-exp(-x)).*(exp(-x)-exp(x-z)).*(1-exp(x-z))./(2-exp(-x)-exp(x-z)).^2 ...
% %     + 1./(2-exp(-x)-exp(x-z)).*((2*exp(-x)-1).*(1-exp(-x))+(2-exp(x)).*(exp(-x)-exp(-z))) );
% % zz = 1:0.001:10;
% % 2 steps
% % Rtheory = (1-exp(-2*z))./(4*(1-exp(-1*z)));
% % Ttheory = z.*exp(-z)./(2*(1-exp(-1*z)));
% % three steps
% % Rtheory=[];Ttheory=[];
% for ii = 1:length(zz)
%     [Rtheory(ii),eb(ii)] = quadgk(@(x)pr3steps(x,zz(ii)),0,zz(ii));
% %     Ttheory(ii) = quadgk(@(x)pt3steps(x,zz(ii)),0,zz(ii));
% end
% figure;
% hold on;
% plot(zz,Rtheory,'b' )
% % % plot(zz,Ttheory,'r' )
% plot(opticalDepths, rCalc, 'b.')
% % plot(opticalDepths, tCalc,'r.')
% % % % grid on
% % ylim([0 1])
% % legend('Probability of reflection (theory)',  'Probability of transmission (theory)',...
% %     'Probability of reflection (Calculated)','Probability of transmission (Calculated)')
% % xlabel('Optical depth')
% % ylabel('R or T')
% % grid on
% % title(sprintf('Theoretical vs calculated probability of reflection/transmission at step 2'))

mask = (results.displacements<0 & results.nSteps == 2);
histogram( results.distances(mask),'normalization', 'pdf' );hold on;
plot(0:0.01:10,gampdf(0:0.01:10,2,1))



