z=1;
rng(645651)
d = struct;s1 = sim.iid.MonteCarloFast(d, 300000, 'Exponential', 1);
s2 = sim.iid.MonteCarloFast(d, 300000, 'Exponential', 1);
p = [];p2=[];
for ii = z
    mask = (s1.result<=ii & s2.result>s1.result);
    p(ii)=nnz(mask)/length(s1.result);
    m2 = s1.result<=ii;
    p2(ii) = nnz(s2.result(m2)>=s1.result(m2))/nnz(m2);
end
zs = 0:0.01:10;
figure;plot(z,p,'b.')
hold on;plot(zs, 0.5*(1-exp(-2*zs)), 'b')
plot(z,p2,'r.')
plot(zs, 0.5*(1-exp(-2*zs))./(1-exp(-zs)), 'r')

[y,x] = ecdf( s1.result(mask)+s2.result(mask) );
c = cdfGuess(x,1);
plot(x,y);hold on;plot(x,c/c(end),'k');

function y=cdfGuess(x,z)
prefactor = (1-exp(-2*x));%/(1-exp(-z));
y = zeros(size(x));
y(x<z) = 1-exp(-x(x<z)).*(x(x<z)+1);
y(x>=z) = (1-exp(-z).*(z+1))+z*(exp(-z)-exp(-x(x>=z)));
y=y.*prefactor;
end