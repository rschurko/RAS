function [RAS, kernel, relax, sys]=Sim_TSVD_R2_RAS_v4(sys,noise,alpha,lam)
% Simulated 2D R2 Relaxation-Assisted Separation (RAS)  

%%% MJJ 2017 %%%
%%% AA  2020 %%%

%-%-%-%-%-%-%-%-%-%-%-% NOTES %-%-%-%-%-%-%-%-%-%-%-% 
%l1 and l2 norm - elastic net (EN)

%-%-%-%-% Experimental Parameters %-%-%-%-%
sys.mg=100;                 %Number of echoes.                INPUT
sys.echosize=512;          %Number of points per echo.       INPUT
sys.dwell=0.001;           %Dwell time (ms).                 INPUT
sys.Tx=0;                  %Set the transmitter offset frequency in kHz

kernel.NumRates=1000;       %Number of R2 relaxation rates to use in the fitting INPUT
kernel.LowerLimit=0.001;      %Lower limit of R2 fit                         INPUT
kernel.UpperLimit=10;      %Upper limit of R2 fit                         INPUT 

relax.T2A = 5; %Simulated T2A in ms
relax.T2B = 30; %Simulated T2B in ms

kernel.e = 3;    %denoising - specify # of singular values to keep - 0 for no denoising
sys.noise = noise;  %Simulated noise

kernel.alpha=alpha;%0.0; %l2 std(real(sys.Spec(1:50,1)));          %Regularization parameter.                     INPUT                    
kernel.lambda=lam;%0.0; %l1 
%kernel.r = 17;    %r number of singular values to keep
%-%-%-%-% %-%-%-%-% %-%-%-%-%

sys.Frequency=transpose(-(1/sys.dwell)/2+sys.Tx:(1/sys.dwell)/sys.echosize:(1/sys.dwell)/2+sys.Tx-(1/sys.dwell)/sys.echosize);

%-%-%-%-% Construct Relaxation Structure %-%-%-%-%
sys.echolength=sys.echosize*sys.dwell;      %Length of one echo in ms.

%Define the taus (t1 points) from the echo tops and store as a column vector (i.e., numel(taus),1)).
%relax.taus=sys.tau;%transpose([sys.echolength/2:sys.echolength:sys.echolength*sys.mg-sys.echolength/2]);  
relax.taus=transpose(sys.echolength/2:sys.echolength:sys.echolength*sys.mg-sys.echolength/2);  
%sys.Acq_Time=sys.mg*sys.echolength; %Total FID time
%Simulated T2 relaxation for sample A and sample B
relax.Mxy_SampleA=transpose(exp(-1*relax.taus/relax.T2A));
relax.Mxy_SampleB=transpose(exp(-1*relax.taus/relax.T2B));

% Calculation for 1D R2 decay curves for the mixture 
relax.Mxy_Mixture=transpose(relax.Mxy_SampleA+relax.Mxy_SampleB);                   % 1D T2 curve for the Mixture

% Calculate the 2D T2 datasets (akin to the datasets present at step 2 of
% the R2 RAS procedure - see Figure 1b in the manuscript).
relax.T2_DataSetA=sys.A*relax.Mxy_SampleA;                        %Full T2 2D dataset for species A.
relax.T2_DataSetB=sys.B*relax.Mxy_SampleB;                        %Full T2 2D dataset for species B
relax.T2_DataSetMix=relax.T2_DataSetA+relax.T2_DataSetB;          %Full T2 2D dataset for the mixture
relax.T2_DataSetMix=relax.T2_DataSetMix/max(relax.T2_DataSetMix,[],'All');

% Adding White Gaussian Noise to the T2 Datasets and 1D Overlapping Lineshapes 
relax.noise_signal_Dataset_sum=sum(sys.noise*randn(size(transpose(relax.T2_DataSetMix)))); 
relax.noise_signal_Dataset=sys.noise*randn(size(relax.T2_DataSetMix));
relax.T2_DataSetMix_noise=((relax.T2_DataSetMix+relax.noise_signal_Dataset)); 
relax.stdn = std(relax.T2_DataSetMix_noise(1:100,1)) ;
relax.snr = max(relax.T2_DataSetMix_noise,[],'all')/relax.stdn;
relax.snrT2 = max(relax.T2_DataSetMix_noise,[],'all')/std(relax.T2_DataSetMix_noise(1,:));
%-%-%-%-% %-%-%-%-% %-%-%-%-%

%-%-%-%-% Construct Kernel Structure %-%-%-%-%
kernel.R=logspace(log10(kernel.LowerLimit),log10(kernel.UpperLimit),kernel.NumRates);   % R2 relaxation rates to be used in kernel. 
%kernel.R=linspace((kernel.LowerLimit),(kernel.UpperLimit),kernel.NumRates);   % R2 relaxation rates to be used in kernel. 
tic
kernel.exponential=-relax.taus*kernel.R;
kernel.M=zeros(sys.mg,kernel.NumRates);
kernel.M=1.0*exp(kernel.exponential);  

%%%Maximum entropy r selection%%%
[m , n] = size(kernel.M);
l = min([m,n]);
s = svd(kernel.M);

for i = 1:l
    ss(i) = s(i)^2;
end
ss = sum(ss);

for i = 1:l
    S(i) = s(i)^2/ss;
    temp(i) = S(i)*log(S(i));
    Etemp(i) = (-1/log(l))*temp(i);
end    
temp = sum(temp);

E = (-1/log(l))*temp;

for i = 1:l
    delE(i) = E - ( (E*log(l) + S(i)*log(S(i)) )/ log(l-1)) ;
end

c = mean(delE);
d = std(delE);

[notr, kernel.r ] = min( delE - (c+d));

%%%Done entropy section - use r for TSVD%%%

[kernel.U, kernel.s, kernel.V]=svd(kernel.M);                                      % Calculate the SVD of the relaxation kernel - U and V are the left
                                                                                    % and right singular vectors and s contains the singular values.
%Mz=[Mz;lambda*eye(min(size(Mz)),max(size(Mz)))];                                   % Tikhonov regularization with L=Identity 

kernel.U = kernel.U(:,1:kernel.r);
kernel.s = kernel.s(1:kernel.r,1:kernel.r);
kernel.V = kernel.V(:,1:kernel.r);

kernel.Mxy = kernel.s*transpose(kernel.V); %Form compressed kernel

kernel.L = get_l(kernel.NumRates,2);                  % Compute the discrete second order derivative.
%kernel.L = eye(kernel.NumRates-1,kernel.NumRates);  %Identity for elastic net reg
kernel.MxyL =[kernel.Mxy;kernel.alpha*kernel.L];   % Tikhonov regularization with L=2nd discrete deriv. op.

% The 2D T2 dataset for the overlapping lineshapes
%relax.T2_DataSetMix_noise=sys.Spec; %((relax.T2_DataSetMix+relax.noise_signal_Dataset)); %%%This is where the experimental data should go
if kernel.e > 0
    relax.T2_DataSetMix_noise = denoise(relax.T2_DataSetMix_noise,kernel.e);
end
relax.Signal=abs(relax.T2_DataSetMix_noise');
relax.stdnpca = std(relax.T2_DataSetMix_noise(1:100,1)) ;
relax.snrpca = max(relax.Signal,[],'all')/relax.stdnpca;
relax.snrpcaT2 = max(relax.Signal,[],'all')/std(relax.T2_DataSetMix_noise(1,:));
% Calculate approximate S/N ratios - NOT USED
%S2N_Dataset=mean(sum(T2_DataSetMix').^2)/mean(noise_signal_Dataset_sum.^2);
%S2N_Mixture=mean(Mixture.^2)/mean(noise_signal_Mixture.^2);
%-%-%-% %-%-%-% 


%-%-%-%-% %-%-%-%-% NNTF routine for each frequency point in the mixture %-%-%-%-% %-%-%-%-%
%relax.T2_DataSetMix_noise=transpose(relax.T2_DataSetMix_noise(1:max(size(sys.Frequency)),1:sys.mg));       % Remove the concatenated regularization section of the data matrix.

% Vertically concatenate the T2 Dataset containing the added noise with 
% zeros so this matrix and the regularization matrix are the appropriate sizes 
% for the non-negative LS MATLAB function (see Eq. 12 in the Theory Section). 
%relax.T2_DataSetMix_noise=[relax.T2_DataSetMix_noise;zeros(size(kernel.L,1),max(size(sys.Frequency)))];  % Tikhonov regularization with L=2nd discrete deriv.

% For each frequency point in the T2 2D dataset, find the R2 relaxation rates
% Initialize the R2 RAS dataset
RAS.spec = zeros(kernel.NumRates,max(size(sys.Frequency)));
dummyspec = zeros(kernel.NumRates,max(size(sys.Frequency)));
%SignalCompress = zeros(2*kernel.r,max(size(sys.Frequency)));
%SignalCompress = zeros(kernel.r,max(size(sys.Frequency)));
%tic

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
%RAS.lam_max = lambda_max;
%-%-%-%-% %-%-%-%-% %-%-%-%-% %-%-%-%-% 

%-%-%-%-% %-%-%-%-% Plotting 2D R2 RAS Spectra %-%-%-%-% %-%-%-%-%
figure(1)
% subplot(2,2,1);
% mesh(sys.Frequency,relax.taus,relax.Signal);
% title('Sample Mixture');
% xlabel('Spectrum Frequency (kHz)');
% %xlim([-120 80]);
% ylabel('Relaxation Time (ms)');
% set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial');
% view([55,15])

subplot(1,2,1);
contour(sys.Frequency,kernel.R,RAS.spec,linspace( 0.0*max(RAS.spec,[],'All'),max(RAS.spec,[],'All'),30));
%contour(sys.Frequency,1./kernel.R,RAS.spec,25);
title('2D RAS Spectrum');
xlabel('Spectrum Frequency (kHz)');
ylabel('R_2 (ms-1)');
%ylim([1,60]);
xlim([-270 180]);
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial','YScale','log');
colormap(jet(25));

subplot(1,2,2);
plot(kernel.R,sum(RAS.spec(:,220:328),2));
%ylim([0 3])
title('R_2 Projection');
xlabel('R_2 (ms-1)');
set(gca, 'XDir','reverse','TickDir','out','FontSize',14,'FontName','Arial','XScale','log');
%subplot(2,2,3);
%plot(sum(RAS.spec(200:400,:),1));
end