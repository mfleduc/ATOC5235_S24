clear variables;close all
x = -25:0.01:25;
N = 10.^(0:2);
f = zeros(length(x), length(N));
legendStr = {};
for ii = 1:length(N)
    f(:,ii) = GenLap(x,N(ii)-1/2 );
    legendStr{ii} = sprintf('N=%d',N(ii));
end
figure;plot(x,f')
grid on
legend(legendStr)
ylabel('f_{S_N}(x)')
xlabel('x')
title('Generalized Laplace PDFs for N=1,10,100')
xlim([min(x),max(x)])
function f = GenLap(x,nu)
%% Generalized Lapplace distribution
numer = abs(x).^nu .* besselk(nu, abs(x));
denom = 2^nu*sqrt(pi)*gamma(nu+0.5);
f = numer/denom ;

end