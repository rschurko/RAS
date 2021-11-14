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

subplot(1,3,1);
%mesh(squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2)));
contour(Freq,R2Rates,squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2))',40);
title('R2 RAS Compression')
xlabel('Frequency (kHz)');
%xlim([-125 75])
ylabel('R2 s^-^1');
set(gca,'XDir','reverse','YScale', 'log');

subplot(1,3,2);
%mesh(RatesDistributionUnStacked(:,:,j));
%mesh(squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3)));
contour(Freq,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3))',40);
title('R1 RAS Compression')
xlabel('Frequency (kHz)');
ylabel('R1 s^-^1');
%xlim([-125 75])
set(gca,'XDir','reverse','YScale', 'log');

subplot(1,3,3);
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),40);
%contour(R2Rates,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(RatesDistributionUnStacked,[],'all'),max(RatesDistributionUnStacked,[],'all'),70));
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(RatesDistributionUnStacked,[],'all'),max(RatesDistributionUnStacked,[],'all'),15));
contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),40));
colormap(jet(25));
title('R1-R2 Map')
xlabel('R2 index');
ylabel('R1 index');

%set(gca,'XScale','log','YScale', 'log');
set(gcf, 'Position',  [10, 80, 700, 600])

%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[0.2*max(RatesDistributionUnStacked,[],'all') max(RatesDistributionUnStacked,[],'all')]);
j = 1;
k = 150;
l = 300;

figure(2)
subplot(1,2,1);
plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,1:k,1:300),2)),2));
title('Pattern 1')
xlabel('Frequency (kHz)');
%xlim([-75 75]);
set(gca,'XDir','reverse');

subplot(1,2,2);
plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,k:l,200:280),2)),2));
title('Pattern 2')
xlabel('Frequency (kHz)');
%xlim([-150 100]);
set(gca,'XDir','reverse');
set(gcf, 'Position',  [710, 80, 600, 600])
