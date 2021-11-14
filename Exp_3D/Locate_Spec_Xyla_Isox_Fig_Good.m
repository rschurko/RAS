MinRateR1 = 0.05;                % 1/s
MaxRateR1 = 50;              % 1/s

MinRateR2 = 1;              % 1/s
MaxRateR2 = 500;              % 1/s

NumRatesR1 = 300;              % Number of R1 rates
NumRatesR2 = 300;             % Number of R2 rates
R1Rates = logspace(log10(MinRateR1),log10(MaxRateR1),NumRatesR1); % R1 Vector
%R1Rates = linspace((MinRateR1),(MaxRateR1),NumRatesR1);
R2Rates = logspace(log10(MinRateR2),log10(MaxRateR2),NumRatesR2); % R2 Vector

load('Xyla_Isox_3D_a1e-4_lam5.mat')
a=103;
Tx = 0;
scale = linspace(-250+Tx,250+Tx,512);
Freq = linspace(scale(198),scale(300),a);


i = 1; %R1 lower bound (y)
i2 = 300; %R1 upper bound (y)

j = 1; %R2 lower bound (x)
j2 = 300; %R2 upper Bound (x)

lvl = 0.00; % percentage of lower threshold in R1/R2 Map

%a = find(R1Rates(i));
%b = find(R2Rates(j));

% subplot(1,3,1);
% %mesh(squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2)));
% contour(Freq,R2Rates,squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2))',40);
% title('R2 RAS Compression')
% xlabel('Frequency (kHz)');
% %xlim([-125 75])
% ylabel('R2 s^-^1');
% set(gca,'XDir','reverse','YScale', 'log');
% 
% subplot(1,3,2);
% %mesh(RatesDistributionUnStacked(:,:,j));
% %mesh(squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3)));
% contour(Freq,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3))',40);
% title('R1 RAS Compression')
% xlabel('Frequency (kHz)');
% ylabel('R1 s^-^1');
% %xlim([-125 75])
% set(gca,'XDir','reverse','YScale', 'log');
% 
% subplot(1,3,3);
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),40);
%contour(R2Rates,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(RatesDistributionUnStacked,[],'all'),max(RatesDistributionUnStacked,[],'all'),70));
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(RatesDistributionUnStacked,[],'all'),max(RatesDistributionUnStacked,[],'all'),15));
contour(R1Rates,R2Rates,squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),40));
colormap(jet(25));
title('R1-R2 Map')
xlabel('R2 (s-1)');
ylabel('R1 (s-1)');
%set(gca,'TickDir','out','FontSize',14,'FontName','Arial','XScale','log','YScale','log','XTick', [1,50,100]);
set(gca,'XScale','log','YScale', 'log','TickDir','out','FontSize',14,'FontName','Arial','XTick',[0.1,0.5,1,5,10,20],'YTick',[0.1,0.5,1,5,10,20,50,100]);

set(gcf, 'Position',  [10, 80, 700, 600])
%set(gca,'TickDir','out','FontSize',14,'FontName','Arial');
%xlim([1.5 30])
%ylim([2 60])

%set(gca,'XScale','log','YScale', 'log');
set(gcf, 'Position',  [10, 80, 700, 600])

%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[0.2*max(RatesDistributionUnStacked,[],'all') max(RatesDistributionUnStacked,[],'all')]);

k = 120;

figure(2)
subplot(1,2,1);
plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,1:k,1:250),2)),2));
title('Pattern 1')
xlabel('Frequency (kHz)');
%xlim([-75 75]);
set(gca,'XDir','reverse');

subplot(1,2,2);
plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,160:190,220:300),2)),2));
title('Pattern 2')
xlabel('Frequency (kHz)');
%xlim([-150 100]);
set(gca,'XDir','reverse');
set(gcf, 'Position',  [710, 80, 600, 600])
