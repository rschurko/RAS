clear;
load('deR1.mat')
sys.Spec = despecR1; %Pre-processed denoised IR data

sys.tau = logspace(log10(0.001),log10(3), 32)'; %Measured IR delays (VDList)
% sys.tau = 0; %%VDlist read syntax

[RAS, kernel, relax, sys] = TSVD_R1_RAS_v4(sys,0.01,0); %Run R1 RAS

Tx = 0;
freq = linspace(-250+Tx,250+Tx,512);

%%Plotting Spectra
%Can use contour(RAS.spec) or mesh(RAS.spec) to find indices to separate:
a = 5;
b = 530;
c = 900;

specA = sum(RAS.spec(a:b,:),1);
specB = sum(RAS.spec(b:c,:),1) ;

figure(2)
subplot(1,2,1)
plot(freq,(specA))
%xlim([-40 30])
xlabel('Frequency (kHz)')
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');

subplot(1,2,2)
plot(freq,specB)
%xlim([-120 120])
xlabel('Frequency (kHz)')
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');