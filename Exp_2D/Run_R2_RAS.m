clear;
load('deR2.mat')
sys.Spec = despecR2; %Pre-processed denoised IR data
load('t2Cl.mat')
sys.tau = t2;  %Pre-calculated echo times

[RAS, kernel, relax, sys] = TSVD_R2_RAS_v4(sys,0.05,1e-2); %Run R2 RAS
%R2 behavior is strange for Isox, l1 norm can help here

Tx = 0;
freq = linspace(-250+Tx,250+Tx,512);

%Can use contour(RAS.spec) or mesh(RAS.spec) to find indices to separate:
a = 1;
b = 815;
c = 950;

specA = sum(RAS.spec(a:b,:),1);
specB = sum(RAS.spec(b:c,:),1);

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