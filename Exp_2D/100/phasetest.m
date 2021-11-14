
dum = sp(:,1);

[ph0, ph1, ph2, off] = autophase(dum,10,0);

out = phi_v3(dum,159.5, 18360.8, -2091.0, -19);

plot(real(out))
hold on
plot(abs(dum))
hold on
plot(abs(out) - real(out) + 1e7);
hold off
%xlim([250 750])
%ylim([0 8*1e7])