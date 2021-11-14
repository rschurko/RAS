function [phase0, phase1, phase2, bestoff] = autophase(spec1,iter,do2)
%Automatic phasing routine for up to 2nd-order phasing.
%do2 is either = 0 (find optimal ph2 and offset) or 1 (do not search ph2 or offset)
%MJJ + ARA 2021

np = numel(spec1);
i = linspace(0,np,(np+1));     % ph1 variables (in pts)
j = linspace(0,360,(500));    % ph0 variables  (in degrees)
l = linspace(400*90,-400*90,(500));  % ph2 variables (in degrees)
ph0_1 = 0;
ph2_1 = zeros(1,iter+1);
BestOffset = zeros(1,iter+1);
offsets = np:-4:-np ;

%% Start Iterations - First Leg
for kk=1:iter
% Find first order phase    
    Iter1{kk} = sum(real(phi_v3(spec1,ph0_1(kk),360*i,ph2_1(kk),round(BestOffset(kk)))));
    M_ph1(kk) = max(abs(Iter1{kk}));
    ph1_1(kk) = i(find(abs(Iter1{kk}) == M_ph1(kk)));
    
    % Find zero order phase
    Iter3{kk} = sum(real(phi_v3(spec1,j,360*ph1_1(kk),ph2_1(kk),round(BestOffset(kk)))));
    M_ph0(kk) = max(real(Iter3{kk}));
    a = find(real(Iter3{kk}) == M_ph0(kk));
    if numel(a) > 1
        a = a(1);
    end
    ph0_1(kk+1) = j(a);
    
    if do2 == 0
        %Find second order phase
        for mm = 1:numel(offsets)
            Iter2{kk,mm} = (sum(real(phi_v3(spec1,ph0_1(kk+1),360*ph1_1(kk),l,(offsets(mm))))));
            M_ph2(kk,mm) = max(abs(Iter2{kk,mm}));
        end
        BestOffsetIndex = find(abs(M_ph2(kk,:)) == max(abs(M_ph2(kk,:))));
        if numel(BestOffsetIndex) > 1
            BestOffsetIndex = BestOffsetIndex(1);
        end
        BestOffset(kk+1) = offsets(BestOffsetIndex);
        ph2_1(kk+1) = l(find((abs(Iter2{kk,BestOffsetIndex})) == abs(M_ph2(kk,round(BestOffsetIndex)))));
    end
    
    area1(kk) = (sum(real(phi_v3(spec1,ph0_1(kk+1),360*ph1_1(kk),ph2_1(kk+1),round(BestOffset(kk+1))))));
    
Highestkk1 = find((area1) == max(area1));

a = Highestkk1(1);    
ph0 = ph0_1(a+1); ph2 = ph2_1(a+1); ph1 = ph1_1(a);
if kk>50
    i = ph1 + linspace(-64/51,64/51,round(500/(51))); % form of dropout...
else
    i = ph1 + linspace(-64/(kk+1),64/(kk+1),round(500/(kk+1)));
end
j = ph0 + linspace(0/(kk+1),360/(kk+1),round(500/(kk+1)));
l = ph2 + linspace(400*90/(kk+1),-400*90/(kk+1),(round(500)));
offsets = BestOffset(a+1)+linspace(128,-128,round(256));
Phases{kk,1} = ph0; Phases{kk,2} = ph1; Phases{kk,3} = ph2; Phases{kk,4} = BestOffset(a+1); 

disp(strcat(num2str(kk),'/ ', num2str(iter)))
end
plot(real(flipud(((phi_v3(spec1,ph0_1(a+1),ph1_1(a)*360,ph2_1(a+1),round(BestOffset(a+1))))))), 'k'); 
hold on
plot(flipud(abs(spec1)), 'r--'); 
hold off
legend('Phased','Magnitude')
phase0 = ph0_1(a+1); phase1 = ph1_1(a)*360; phase2 = ph2_1(a+1);
bestoff = round(BestOffset(a+1));

fprintf('ph0, ph1, ph2, off = %.1f, %.1f, %.1f, %.1d\n',phase0,phase1,phase2,bestoff)
end

