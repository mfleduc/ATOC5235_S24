%%  run file to test convergence of the Monte Carlo algorithm
clear variables;close all;clc;
rng(645651)
possiblePaths = {'/homes/male7736/Desktop/Research/Monte-Carlo-Code',...
    'C:\Users\Test\Documents\Research\Monte-Carlo-Code'} ;
addpath( possiblePaths{ispc+1} );%Monte Carlo Code
input=struct();
%% Input parameters
nPhotons = ceil(10.^[1:0.1:6]); %Number of photons to simulate
meanFreePath = 1; %Mean free path: Mean distance traveled before interaction 
%% We have to see how long it takes to converge to an exponential distribution
%Test convergence by taking the sup norm of the difference between the
%empirical CDF and the real distribution function
yecdfs = {};
xecdfs = {};
trueDist = makedist('Exponential','mu',meanFreePath);
truey = {} ;
ksTestErr = [];
ksTestConv = [] ;
for kk = 1:length(nPhotons)
    out = sim.iid.MonteCarloFast( struct.empty(), nPhotons(kk), 'Exponential',meanFreePath ) ;
    [yecdfs{kk},xecdfs{kk}] = ecdf( out.result );
    truey{kk} = trueDist.cdf(xecdfs{kk});
    ksTestErr(kk) = max(abs(truey{kk}-yecdfs{kk}));
end

figure;loglog(nPhotons, ksTestErr)
grid on;
xlabel('Number of photons')
ylabel('||F_{est}-F_{true}||_{\infty}')
title('Convergence in the sup-norm of Monte-Carlo method')

ndx = find( ksTestErr < 1e-2  );

figure;
plot( xecdfs{ndx(1)}, yecdfs{ndx(1)}, 'b' )
hold on
plot(0:0.01:max(xecdfs{ndx(1)}), trueDist.cdf(0:0.01:max(xecdfs{ndx(1)})),'r')
grid on
legend('F_{est}','F_{true}')
xlim([0 5])
xlabel('x')
ylabel('F(x)')
title('Estimated vs true CDF, n \approx 8000')