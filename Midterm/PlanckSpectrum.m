function spectrum = PlanckSpectrum(T,L)
c=3e8;
h=6.626e-34;
kb = 1.381e-23;
spectrum = zeros(size(L));

specDenom = L.^5.*(exp(h*c./(kb*L*T))-1);
specNum = 2*h*c^2;
spectrum = specNum./specDenom;
end