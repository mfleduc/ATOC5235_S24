function result = ThreeDSimulate(input)
%% result = simulate(input)
%   Code relying on the github repository for Monte Carlo code located at 
%   https://github.com/mfleduc/Monte-Carlo-Code for solving the radiative 
%   transfer equations. Before starting, this repository should be added to
%   the MATLAB path.
%   This code solves the radiative transfer equations to determine the
%   reflection, transmission, and absorption coefficients in a 3D layer
%   with infinite horizontal extent
nPhotonsRemoved = 0;
maxSteps = max(10, ceil(1e2*input.opticalDepth));
distFn = struct();
g = input.asymmetry;
switch input.scatteringProbs
    case 'Henyey-Greenstein'
        distType = 'pdf';
        
        distFn.x = 0:0.01:pi;
        distFn.y = 1/(2)*(1-g^2)./(1+g^2-2*g*cos(distFn.x)).^(3/2);
    otherwise
        distType = input.scatteringProbs.type;
        distFn.x = input.scatteringProbs.x;
        distFn.y = input.scatteringProbs.y ;
end
% distFn.x = 0:1/(100*input.mfp):50*input.mfp;
% distFn.y = exp( -1*input.mfp*distFn.x );
R=0;T=0;A=0; %Reflected, transmitted, absorbed photon count
bins = linspace(0,input.opticalDepth, 1e3);
binCounter = zeros(2,length(bins));%Down,up
basisVectors = eye(3);
for ph = 1:input.Nphotons
    taus = zeros(3, maxSteps+1);
    stepSizes = sim.iid.MonteCarloFast( distFn,...
        maxSteps+1, 'Exponential',input.mfp );
    taus(:,1) = [0;0;stepSizes.result(1)];
    tauT = taus(3,1);
    %Simulate all the steps right at the start to prevent wasting
    %time making repeated calls to sim.iid.MonteCarloFast()
    phi = 2*pi*rand( 1, maxSteps+1 );%Will convert these into angles
    %in the while loop
    rndForAbsorption = rand(1, maxSteps+1);
    thetas = sim.iid.MonteCarloFast( distFn,...
        maxSteps+1,distType );%Storing the 
    %angles
    ct = cos(thetas.result);
    st = sin(thetas.result);
    cp = cos(phi);
    sp = sin(phi);
    absFlag = 0;
    cnt=0;
    
    while ~absFlag&&tauT>0&&tauT<input.opticalDepth &&...
            cnt<maxSteps
        % While the photon is not absorbed, the photon is in the
        % cloud, and we don't have a ton of bouncing around (this
        % is to avoid an accidental infinite loop)
        cnt=cnt+1;
        if rndForAbsorption(cnt)<=input.ssa
            Omega = basisVectors*[st(cnt)*cp(cnt),st(cnt)*sp(cnt),ct(cnt)]';
            taus(:,cnt+1) = taus(:,cnt)+Omega*stepSizes.result(cnt+1) ;
            M = [[ct(cnt)*cp(cnt), -1*sp(cnt), st(cnt)*cp(cnt)];...
                 [ct(cnt)*sp(cnt), cp(cnt),    st(cnt)*sp(cnt)];...
                 [-1*st(cnt),      0,           ct(cnt)       ]];
            basisVectors = basisVectors*M ;
            tauT = taus(3,cnt+1);
        else
            absFlag=1;
        end
    end
    if absFlag
        A=A+1;
    elseif tauT<0
        R=R+1;
    elseif tauT>=input.opticalDepth
        T=T+1;
    else
        nPhotonsRemoved = nPhotonsRemoved+1;
    end
end

R = R/(input.Nphotons-nPhotonsRemoved);
T = T/(input.Nphotons-nPhotonsRemoved);
A = A/(input.Nphotons-nPhotonsRemoved);
result = struct();
result.errCode = 0;
result.RTA = [R,T,A];

end

