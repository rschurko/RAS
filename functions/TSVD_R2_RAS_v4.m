function [RAS, kernel, relax, sys]=TSVD_R2_RAS_v4(sys,alpha,lam)

%%% MJJ 2017 %%%
%%% AA  2020 %%%

%-%-%-%-%-%-%-%-%-%-%-% NOTES %-%-%-%-%-%-%-%-%-%-%-% 
% EXP Sparse NNTF w/ Tikhonov (l2) and l1 penalties
% selects TSVD r parameter with maximum entropy criterion

%-%-%-%-% Experimental Parameters %-%-%-%-%
sys.Spec = real(sys.Spec);
sys.Spec = sys.Spec / max(sys.Spec,[],'All');
[sys.echosize,sys.mg] = size(sys.Spec);

sys.dwell=0.002;           %Dwell time (ms).                 INPUT
sys.Tx=0;                  %Set the transmitter offset frequency in kHz

kernel.NumRates=1000;       %Number of R2 relaxation rates to use in the fitting INPUT
kernel.LowerLimit=0.1;      %Lower limit of R2 fit                         INPUT
kernel.UpperLimit=1000;       %Upper limit of R2 fit                         INPUT 

kernel.alpha = alpha; %l2 std(real(sys.Spec(1:50,1)));          %Regularization parameter.                     INPUT                    
kernel.lambda = lam; %l1                   INPUT                    
%-%-%-%-% %-%-%-%-% %-%-%-%-%

sys.Frequency=transpose(-(1/sys.dwell)/2+sys.Tx:(1/sys.dwell)/sys.echosize:(1/sys.dwell)/2+sys.Tx-(1/sys.dwell)/sys.echosize);

%-%-%-%-% Construct Relaxation Structure %-%-%-%-%
sys.echolength=sys.echosize*sys.dwell;      %Length of one echo in ms.

%sys.tau is the times of the echo tops:
relax.taus = sys.tau';
%-%-%-%-% %-%-%-%-% %-%-%-%-%

%-%-%-%-% Construct Kernel Structure %-%-%-%-%
kernel.R=logspace(log10(kernel.LowerLimit),log10(kernel.UpperLimit),kernel.NumRates);   % R2 relaxation rates to be used in kernel. 
%kernel.R=linspace((kernel.LowerLimit),(kernel.UpperLimit),kernel.NumRates);   % R2 relaxation rates to be used in kernel. 

kernel.exponential=-relax.taus*kernel.R;
kernel.M=zeros(sys.mg,kernel.NumRates);
kernel.M=1*exp(kernel.exponential);      %T2 exponential decay for kernel    

%%%Maximum entropy r selection%%%
[m , n] = size(kernel.M);
l = min([m,n]); s = svd(kernel.M);

for i = 1:l
    ss(i) = s(i)^2;
end
ss = sum(ss);

for i = 1:l
    S(i) = s(i)^2/ss;
    temp(i) = S(i)*log(S(i));
    Etemp(i) = (-1/log(l))*temp(i);
end    
temp = sum(temp); E = (-1/log(l))*temp;

for i = 1:l
    delE(i) = E - ( (E*log(l) + S(i)*log(S(i)) )/ log(l-1)) ;
end
c = mean(delE);
d = std(delE);

[notr, kernel.r ] = min( delE - (c+d));
%%%Done entropy section - use r for TSVD%%%
[kernel.U, kernel.s, kernel.V]=svd(kernel.M);  % Calculate the SVD of the kernel                                                                                % and right singular vectors and s contains the singular values.

kernel.U = kernel.U(:,1:kernel.r);
kernel.s = kernel.s(1:kernel.r,1:kernel.r);
kernel.V = kernel.V(:,1:kernel.r);

kernel.Mxy = kernel.s*kernel.V'; %Form compressed kernel

kernel.L = get_l(kernel.NumRates,2);                  % Compute the discrete second order derivative.
%kernel.L = eye(kernel.NumRates-1,kernel.NumRates);  %Identity for elastic net reg
kernel.MxyL =[kernel.Mxy;kernel.alpha*kernel.L];   % Tikhonov regularization with L=2nd discrete deriv. op.

% The 2D T2 dataset for the overlapping lineshapes
%relax.T2_DataSetMix_noise=sys.Spec; %((relax.T2_DataSetMix+relax.noise_signal_Dataset)); %%%This is where the experimental data should go
relax.stdn = std(sys.Spec(1:100,1));
relax.Signal=real(sys.Spec');
%-%-%-% %-%-%-% 


%-%-%-%-% %-%-%-%-% NNTF routine for each frequency point in the mixture %-%-%-%-% %-%-%-%-%
% Vertically concatenate the T2 Dataset containing the added noise with 
% zeros so this matrix and the regularization matrix are the appropriate sizes 
% for the non-negative LS MATLAB function (see the Theory Section in paper). 

% For each frequency point in the T2 2D dataset, find the R2 relaxation rates
% Initialize the R2 RAS dataset
RAS.spec=zeros(kernel.NumRates,max(size(sys.Frequency)));

tic
if kernel.lambda == 0
    parfor i = 1:max(size(sys.Frequency))
    SignalCompress(:,i) = [kernel.U'*relax.Signal(:,i); zeros(size(kernel.L,1),1)]
    dummyspec(:,i)=lsqnonneg(kernel.MxyL,SignalCompress(:,i)); 
    disp(['RAS step',num2str(i)]);
    end
else
    parfor i = 1:max(size(sys.Frequency))
    SignalCompress(:,i) = [kernel.U'*relax.Signal(:,i); zeros(size(kernel.L,1),1)]
    %lambda_max(i) = find_lambdamax_l1_ls_nonneg(kernel.MxyL',SignalCompress(:,i));
    dummyspec(:,i)=l1_ls_nonneg(kernel.MxyL,SignalCompress(:,i),kernel.lambda,1e-3,true);
    disp(['RAS step',num2str(i)]);
    end
end
RAS.time = toc;
disp('Finished!');
relax.SignalCompress = SignalCompress;
RAS.spec = dummyspec;
%RAS.resn = RAS.resn/max(RAS.resn);
%-%-%-%-% %-%-%-%-% %-%-%-%-% %-%-%-%-% 

%-%-%-%-% %-%-%-%-% Plotting 2D R2 RAS Spectra %-%-%-%-% %-%-%-%-%
figure(1)

title('Sample Mixture');
xlabel('Spectrum Frequency (kHz)');
ylabel('Relaxation Time (ms)');
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');
view([55,15])

subplot(1,2,1);
contour(sys.Frequency,kernel.R,RAS.spec,linspace( 0.0*max(RAS.spec,[],'All'),max(RAS.spec,[],'All'),100));
title('2D RAS Spectrum');
xlabel('Spectrum Frequency (kHz)');
ylabel('Rate');
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial','YScale','log');
colormap(jet(25));

subplot(1,2,2);
plot(kernel.R,sum(RAS.spec(:,:),2));
xlabel('Rate');
set(gca,'TickDir','out','FontSize',14,'FontName','Arial','XScale','log');
end