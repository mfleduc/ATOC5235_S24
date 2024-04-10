function result = TwoStreamSimulate(input)
%% result = simulate(input)
%   Code relying on the github repository for Monte Carlo code located at 
%   https://github.com/mfleduc/Monte-Carlo-Code for solving the radiative 
%   transfer equations. Before starting, this repository should be added to
%   the MATLAB path.
%  
nPhotonsRemoved = 0;
maxSteps = max(10, ceil(1e2*input.opticalDepth));
distFn = struct.empty();
% distFn.x = 0:1/(100*input.mfp):50*input.mfp;
% distFn.y = exp( -1*input.mfp*distFn.x );
R=0;T=0;A=0; %Reflected, transmitted, absorbed photon count
bins = linspace(0,input.opticalDepth, 1e3);
binCounter = zeros(2,length(bins));%Down,up
travelDists = zeros(1,input.Nphotons);
travelDisps = zeros(1,input.Nphotons);
nSteps = zeros(1,input.Nphotons);
switch input.flags.model
    case 'twostream'
        
        for ph = 1:input.Nphotons
            taus = zeros(1, maxSteps+1);
            stepSizes = sim.iid.MonteCarloFast( distFn,...
                maxSteps+1,'Exponential',input.mfp );
            taus(1) = stepSizes.result(1);
            %Simulate all the steps right at the start to prevent wasting
            %time making repeated calls to sim.iid.MonteCarloFast()
            tmp = rand( 1, maxSteps+1 );%Will convert these into angles
            %in the while loop
            rndForAbsorption = rand(1, maxSteps+1);
            costhetas = ones(1, maxSteps+1);%Storing the cosine of the 
            %angles
            tauT = taus(1);
            binCounter(1,bins<=tauT) = 1+binCounter(1,bins<=tauT);
            absFlag = 0;
            cnt=0;
            lastDirFlag = 0;%0=down,1 = up;
            while ~absFlag&&tauT>0&&tauT<input.opticalDepth &&...
                    cnt<maxSteps
                % While the photon is not absorbed, the photon is in the
                % cloud, and we don't have a ton of bouncing around (this
                % is to avoid an accidental infinite loop)
                cnt=cnt+1;
                if rndForAbsorption(cnt)<=input.ssa
                    taus(cnt+1) = stepSizes.result(cnt+1);
                    if input.flags.isotropicMedium
                        costhetas(cnt+1) = 1-(tmp(cnt+1)>= ...
                            input.scatteringProbs(1))*2;
                    end
                    delTau = taus(cnt+1)*prod( costhetas );
                    ndcs = (bins<max(tauT, tauT+delTau) & ...
                        bins>min(tauT, tauT+delTau));
                    if sign(costhetas(cnt+1))==1
                        %Continuing in the same direction
                        binCounter(lastDirFlag+1, ndcs) =...
                            binCounter(lastDirFlag+1, ndcs)+1;
                    else %Changed direction
                        lastDirFlag = ~lastDirFlag;
                        binCounter(lastDirFlag+1, ndcs) =...
                            binCounter(lastDirFlag+1, ndcs)+1;
                    end
                    
                    tauT = tauT+delTau;
                else
                    absFlag=1;
                end
            end
            if absFlag
                A=A+1;
            elseif tauT<0
                R=R+1;
                nSteps(ph) = cnt+1;
                travelDists(ph) = sum(stepSizes.result(1:cnt+1));
                travelDisps(ph) = tauT ;
            elseif tauT>=input.opticalDepth
                T=T+1;
                nSteps(ph) = cnt+1;
                travelDists(ph) = sum(stepSizes.result(1:cnt+1));
                travelDisps(ph) = tauT ;

            else
                nPhotonsRemoved = nPhotonsRemoved+1;
            end
        end
        
    otherwise
        error('Only two stream model is implemented');
end
FDwnUp = binCounter/(input.Nphotons);
R = R/(input.Nphotons-nPhotonsRemoved);
T = T/(input.Nphotons-nPhotonsRemoved);
A = A/(input.Nphotons-nPhotonsRemoved);
result = struct();
result.errCode = 0;
result.RTA = [R,T,A];
result.FDwnUp = FDwnUp ;
result.distances = travelDists;
result.displacements = travelDisps; 
result.nSteps = nSteps;
end

