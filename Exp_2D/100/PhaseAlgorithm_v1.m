% % Set initial phase parameters 
% clearvars -except PhaseData pieceForPhase
% i = linspace(0,1024,(1001));     % ph1 variables (in pts)
% j = linspace(0,360,(500));    % ph0 variables  (in degrees)
% l = linspace(400*90,-400*90,(500));  % ph2 variables (in degrees)
% ph0_1 = 0;
% ph0_2 = 0;
% ph2_1 = 0;
% ph2_2 = 0;
% BestOffset = 0;
% offsets = [1024:-4:-1024];
% % Set iteration counter
% NumOfIterations = 10;
np = 512;
 for counter=1
%counter=21;
    clearvars -except PhaseData pieceForPhase counter specForMC np spec1
i = linspace(0,np,(np+1));     % ph1 variables (in pts)
j = linspace(0,360,(500));    % ph0 variables  (in degrees)
l = linspace(400*90,-400*90,(500));  % ph2 variables (in degrees)
ph0_1 = 0;
ph0_2 = 0;
ph2_1 = 0;
ph2_2 = 0;
BestOffset = 0;
offsets = [np:-4:-np];
% Set iteration counter
NumOfIterations = 30;
%q = 7;
%counter = 7;
%spec1 = pieceForPhase{1,counter};
%counter = 5; 
%spec1 = spec{14,1};
%spec1 = pieceForPhase{1,counter};
%spec1 = Scheme1_spinHalf_WCPMG_coaddedEchoes; % for scheme1 generation

%% Start Iterations - First Leg
for kk=1:NumOfIterations
        % Real-time for first-order phase
%         for mm=1:numel(i)
%             subplot(1,2,1)
%             plot(real(phi_v2(spec1,ph0_1(kk),360*i(mm),ph2_1(kk))));
%             subplot(1,2,2)
%             plot(mm,sum(real(phi_v2(spec1,ph0_1(kk),360*i(mm),ph2_1(kk)))),'or'); hold on;
%             drawnow;
%             title(num2str(mm))
%         end
    
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
    
    area1(kk) = (sum(real(phi_v3(spec1,ph0_1(kk+1),360*ph1_1(kk),ph2_1(kk+1),round(BestOffset(kk+1))))));
    
Highestkk1 = find((area1) == max(area1));

a = Highestkk1(1);    
ph0 = ph0_1(a+1);
ph2 = ph2_1(a+1);
ph1 = ph1_1(a);
%i = ph1 + [(-20/kk):(40)/(500):(20/kk)];
%j = ph0 + [-(0/(kk^2)):(360/500):360/(kk^2)];
%l = ph2 + [200*90/kk:-22.25/kk:-200*90/kk];
if kk>50
i = ph1 + linspace(-64/51,64/51,round(500/(51))); % form of dropout...
else
    i = ph1 + linspace(-64/(kk+1),64/(kk+1),round(500/(kk+1)));
end
j = ph0 + linspace(0/(kk+1),360/(kk+1),round(500/(kk+1)));
l = ph2 + linspace(400*90/(kk+1),-400*90/(kk+1),(round(500)));
offsets = BestOffset(a+1)+linspace(128,-128,round(256));
Phases{kk,1} = ph0; Phases{kk,2} = ph1; Phases{kk,3} = ph2; Phases{kk,4} = BestOffset(a+1); 
clearvars -except PhaseData specForMC pieceForPhase BestOffset counter offsets M_ph1 M_ph2 M_ph0 spec1 i j l a ph0_1 ph2_1 ph1_1 NumOfIterations Iter1 Iter2 Iter3 area1 Phases kk q

kk
end
PhaseData{1,counter+1} = ph0_1;
PhaseData{2,counter+1} = ph1_1;
PhaseData{3,counter+1} = ph2_1;
PhaseData{4,counter+1} = BestOffset;
PhaseData{5,counter+1} = a;
PhaseData{6,counter+1} = (real(flipud(((phi_v3(spec1,ph0_1(a+1),ph1_1(a)*360,ph2_1(a+1),round(BestOffset(a+1))))))));
end
%figure()
plot(real(flipud(((phi_v3(spec1,ph0_1(a+1),ph1_1(a)*360,ph2_1(a+1),round(BestOffset(a+1))))))), 'k'); 
phase0 = ph0_1(a+1);
phase1 = ph1_1(a)*360;
phase2 = ph2_1(a+1);
bestoff = round(BestOffset(a+1));
% % Find Largest Area
% Highestkk1 = find(area1 == max(area1));
% a = Highestkk1(1);
% ph0 = ph0_1;
% ph1 = ph1_1;
% ph2 = ph2_1;
% %plot(real(flipud(((phi_v3(spec1,ph0(a+1),ph1(a)*360,ph2(a+1),BestOffset(a)))))), 'k'); 
% 
% 
% % ph0(a+1);
% % ph1(a)*360;
% % ph2(a+1);
% % BestOffset(a);
% 
% i = ph1(a) + [-20:40/(500):20 - 40/(500)];
% j = ph0(a+1) + [-50:100/(500):50 - 100/(500)];
% l = ph2(a+1) + [-5*360:10*90/(500):5*360 - 10*90/(500)];
% ph0_1 = ph0(a+1);
% ph2_1 = ph2(a+1);
% BestOffset = BestOffset(a+1);
% clearvars -except BestOffset spec1 i j l ph0_1 ph2_1 NumOfIterations

% %% Start Iterations - Second Leg
% for kk=1:NumOfIterations
%         % Real-time for first-order phase
% %         for mm=1:numel(i)
% %             subplot(1,2,1)
% %             plot(real(phi_v2(spec1,ph0_1(kk),360*i(mm),ph2_1(kk))));
% %             subplot(1,2,2)
% %             plot(mm,sum(real(phi_v2(spec1,ph0_1(kk),360*i(mm),ph2_1(kk)))),'or'); hold on;
% %             drawnow;
% %             title(num2str(mm))
% %         end
%     
% % Find first order phase    
%     Iter1{kk} = sum(real(phi_v3(spec1,ph0_1(kk),360*i,ph2_1(kk),BestOffset(kk))));
%     M_ph1(kk) = max(abs(Iter1{kk}));
%     ph1_1(kk) = i(find(abs(Iter1{kk}) == M_ph1(kk)));
%     
%     % Find zero order phase
%     Iter3{kk} = sum(real(phi_v3(spec1,j,360*ph1_1(kk),ph2_1(kk),BestOffset(kk))));
%     M_ph0(kk) = max(real(Iter3{kk}));
%     a = find(real(Iter3{kk}) == M_ph0(kk));
%     if numel(a) > 1
%         a = a(1);
%     end
%     ph0_1(kk+1) = j(a);
%     
%     %Find second order phase
%     offsets = [120:-1:-120];
%     for mm = 1:numel(offsets)
%         Iter2{kk,mm} = (sum(real(phi_v3(spec1,ph0_1(kk+1),360*ph1_1(kk),l,offsets(mm)))));
%         M_ph2(kk,mm) = max(abs(Iter2{kk,mm}));
%     end
%     BestOffsetIndex = find(abs(M_ph2(kk,:)) == max(abs(M_ph2(kk,:))));
%     if numel(BestOffsetIndex) > 1
%         BestOffsetIndex = BestOffsetIndex(1);
%     end
%     BestOffset(kk+1) = offsets(BestOffsetIndex);
%     ph2_1(kk+1) = l(find((abs(Iter2{kk,BestOffsetIndex})) == M_ph2(kk,BestOffsetIndex)));
%     
%     area1(kk) = sum(real(phi_v3(spec1,ph0_1(kk+1),360*ph1_1(kk),ph2_1(kk+1),BestOffset(kk))));
% kk
% 
% end
% 
% % Find Largest Area
% Highestkk1 = find(area1 == max(area1));
% a = Highestkk1(1);
% ph0 = ph0_1;
% ph1 = ph1_1;
% ph2 = ph2_1;
% %plot(real(flipud(((phi_v3(spec1,ph0(a+1),ph1(a)*360,ph2(a+1),BestOffset(a+1)))))), 'k'); 
% 
% % Find Largest Area
% Highestkk1 = find(area1 == max(area1));
% a = Highestkk1(1);
% ph0 = ph0_1;
% ph1 = ph1_1;
% ph2 = ph2_1;
