%%Plot the L-curve for various alphas
clear;
load('35Cl_SIMPSON_complex.mat')

%alpha = [0.00001, 0.0001, 0.001, 0.01, 0.1];
noise = 0.001; %Amount of simulated noise
alpha = logspace(log10(1e-6),log10(50),11);
lam = 0;

for i =1:numel(alpha)
    [RAS, kernel, relax, sys] = Sim_TSVD_R2_RAS_v4(sys,noise,alpha(i),lam);
    close
for j = 1:512
    res(i,j) = norm( kernel.MxyL*RAS.spec(:,j) - relax.SignalCompress(:,j) );
    diff = (kernel.MxyL*RAS.spec(:,j) - relax.SignalCompress(:,j));
    MSE(i,j) = mean(diff.^2);
    L(i,j) = norm( kernel.L*RAS.spec(:,j) );
    %res(i,j) = norm( kernel.M*RAS.spec(:,j) - relax.Signal(:,j) );
    %L(i,j) = norm( kernel.L*RAS.spec(:,j) );
end
end

%%Sum over pattern frequencies rather than having multiple curves
res = sum(res(:,185:296),2);
L = sum(L(:,185:296),2);
MSE = sum(MSE(:,185:296),2);

plot(log(res),log(L),'.-')
set(gca,'TickDir','out','FontSize',14,'FontName','Arial','XScale','log');%,'YScale','log');
ylabel('log||L*f||')
xlabel('log||K*f - s||')

