%%Run simulated R2-RAS function
%%Open the R2_RAS function to adjust noise and rates
clear;
load('35Cl_SIMPSON_complex.mat')
noise = 0.001; %amount of noise (b/t 0 and 1)
[RAS, kernel, relax, sys]=Sim_TSVD_R2_RAS_v4(sys,noise,0.01,0) ;

%%Select relaxation indices; use mesh(RAS.spec) to find
a = 50;  %1st spectrum start index
b = 480;  %Separatino index
c = 950;  %2nd spectrum end index

%%Plotting
specA = sum(RAS.spec(a:b,:),1);% + sum(RAS.spec(end,:),1);
specB = sum(RAS.spec(b:c,:),1);

subplot(1,2,1)
plot(sys.Frequency,(specA))
xlim([-350 250])
xlabel('Frequency (kHz)')
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');

subplot(1,2,2)
plot(sys.Frequency,specB)
xlim([-350 250])
xlabel('Frequency (kHz)')
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');