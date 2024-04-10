z=1:2;
d = struct;s1 = sim.iid.MonteCarloFast(d, 1e5, 'Exponential', 1);
s2 = sim.iid.MonteCarloFast(d, 1e5, 'Exponential', 1);
p = [];p2=[];
for ii = z
    mask = (s1.result<=ii & s2.result>s1.result);
    p(ii)=nnz(mask)/1e5;
    m2 = s1.result<=ii;
    p2(ii) = nnz(s2.result(m2)>=s1.result(m2))/nnz(m2);
end
zs = 0:0.01:10;
figure;plot(z,p,'b.')
hold on;plot(zs, 0.5*(1-exp(-2*zs)), 'b')
plot(z,p2,'r.')
plot(zs, 0.5*(1-exp(-2*zs))./(1-exp(-zs)), 'r')