%Proc_IR_RAS
Tx = 0;
scale = linspace(-250+Tx,250+Tx,512);
sp1D = sum(rawspecR2(:,1:35),2);

plot(scale, sp1D);
set(gca,'TickDir','out','FontSize',14,'FontName','Arial','XDir','reverse');
xlim([-80 80]);