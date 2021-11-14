%%ARA S-curve code using simulated R2 RAS
clear;
load('35Cl_SIMPSON_complex.mat')

noise = 0.001; %Amount of simulated noise
alpha = 0.001;
lam = logspace(log10(5e-5),log10(1000),10);

for i =1:numel(lam)
    [RAS, kernel, relax, sys]=Sim_TSVD_R2_RAS_v4(sys,noise,alpha,lam(i));
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

res = sum(res(:,185:296),2);
L = sum(L(:,185:296),2);
MSE = sum(MSE(:,185:296),2);

plot((lam),log(res),'.-')
set(gca,'TickDir','out','FontSize',14,'FontName','Arial','XScale','log');%,'YScale','log');
xlabel('\lambda')
ylabel('log||K*f - s||')
