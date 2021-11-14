%%Locate the spectra from the rate distribution object

lvl = 0.00; % percentage of lower threshold in R1/R2 Map
fig2 = 1; %turn on pattern locater with non-zero value


i = 1; %R1 lower bound (y)
i2 = 100; %R1 upper bound (y)
j = 1; %R2 lower bound (x)
j2 = 200; %R2 upper Bound (x)

subplot(1,3,1);
%mesh(squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2)));
contour(Freq,R2Rates,squeeze(sum(RatesDistributionUnStacked(:,i:i2,:),2))',40);
title('R2 RAS Compression')
xlabel('Frequency (kHz)');
%xlim([-125 75])
ylabel('R2 s^-^1');
set(gca,'XDir','reverse','YScale', 'log','TickDir','out');

subplot(1,3,2);
%mesh(RatesDistributionUnStacked(:,:,j));
%mesh(squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3)));
contour(Freq,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,j:j2),3))',40);
title('R1 RAS Compression')
xlabel('Frequency (kHz)');
ylabel('R1 s^-^1');
%xlim([-125 75])
set(gca,'XDir','reverse','YScale', 'log','TickDir','out');

subplot(1,3,3);
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),40);
%contour(R2Rates,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(RatesDistributionUnStacked,[],'all'),max(RatesDistributionUnStacked,[],'all'),70));
contour(R2Rates,R1Rates,squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),40));
%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),linspace(lvl*max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),max(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[],'all'),40));
colormap(jet(25));
title('R1-R2 Map')
%xlabel('R2 index');
%ylabel('R1 index');
xlabel('R2 (s-1)');
ylabel('R1 (s-1)');
set(gca,'XScale','log','YScale', 'log','TickDir','out','FontSize',14,'FontName','Arial');
set(gcf, 'Position',  [10, 80, 700, 600])

%contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)),[0.2*max(RatesDistributionUnStacked,[],'all') max(RatesDistributionUnStacked,[],'all')]);

%%Pattern sep
%%Use mesh() or contour(squeeze(sum(RatesDistributionUnStacked(:,:,:),1)))
%Draw a 'box' around first peak
R21 = 50;
R22 = 115;

R11 = 30;
R12 = 65;

%Draw a 'box' around second peak
R23 = 110;
R24 = 150;

R13 = 65;
R14 = 90;

if fig2 ~= 0
    figure(2)
    subplot(1,2,1);
    plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,R11:R12,R21:R22),2)),2));
    title('Pattern 1')
    xlabel('Frequency (kHz)');
    %xlim([-75 75]);
    set(gca,'XDir','reverse');

    subplot(1,2,2);
    plot(Freq,sum(squeeze(sum(RatesDistributionUnStacked(:,R13:R14,R23:R24),2)),2));
    title('Pattern 2')
    xlabel('Frequency (kHz)');
    %xlim([-150 100]);
    set(gca,'XDir','reverse');
    set(gcf, 'Position',  [710, 80, 600, 600])
end