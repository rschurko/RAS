% R1/R2 Correlation Map - Version 5: 
% 3D Simulated RAS Processing w/ TSVD + EN  + denoising
% ARA 2020
clear;
load('35Cl_SIMPSON_complex.mat')
% Generate Relaxation Data
%t1 = [0:200:200*(31)]*1e-6;      % Number of t1 points sampling T1
t1 = logspace(log10(0.1),log10(3),32);      % Number of t1 points sampling T1
%t2 = [0:500:500*(63)]*1e-6;      % Number of t2 points sampling T2
t2 = linspace(250,250*64,64)*1e-6;      % Number of t2 points sampling T2

R1A = 1.5;                         % R1 Rate Constant, 1/s
R1B = 4;                         % R1 Rate Constant, 1/s
R2A = 300;                         % R2 Rate Constant 1, 1/s  %other ones were 200
R2B = 550;                        % R2 Rate Constant 2, 1/s   %other ones were ca. 600

Data_t1A = 1-2*exp(-t1.*R1A);         % T1 Recovery 
Data_t1B = 1-2*exp(-t1.*R1B);         % T1 Recovery 
Data_t2A = exp(-t2.*R2A);    % T2 Decay
Data_t2B = exp(-t2.*R2B);    % T2 Decay

alpha = 0.0; %0.015;
lambda = 0.0; %0.6;
Noise = 0.0;  %0.014
e = 0;      %denoising; operates w/ t2 dimension
f = 2; %1 or 2 for denoise dimension

%See adjustable kernel params at line 90

% Generate Spectral Data
SiteA = sys.A(137:321); %truncate baseline
SiteB = sys.B(137:321); 

SpectralDataRaw = zeros(numel(SiteA),numel(t1),numel(t2));
T2mix = SiteA*Data_t2A + SiteB*Data_t2B;
T2mix = T2mix/max(T2mix);

for i=1:numel(SiteA)
for   j=1:numel(t1)
SpectralDataRaw(i,j,:) = SiteA(i)*Data_t1A(j).*Data_t2A+SiteB(i)*Data_t1B(j).*Data_t2B;
end
end

SpectralDataRaw = SpectralDataRaw/max(SpectralDataRaw,[],'All');
if Noise == 0
    Noise = 1e-5;
end
SpectralDataRaw = SpectralDataRaw + Noise*randn(numel(SiteA),32,64);
for i=1:numel(t1)
for   j=1:numel(t2)
RawforNoise(:,i,j) = phi_v3(SpectralDataRaw(:,i,j),330.46-0,16606,-21.042,101.7686);
end
end

if e~=0
    if f == 1
        for i=1:numel(t1)
            SpectralDataRaw(:,i,:) = denoise(squeeze(SpectralDataRaw(:,i,:)),e);
        end
    else
        for i=1:numel(t2)
            SpectralDataRaw(:,:,i) = denoise(SpectralDataRaw(:,:,i),e);
        end
    end
end

for i=1:numel(t1)
for   j=1:numel(t2)
SpectralDataRaw(:,i,j) = phi_v3(SpectralDataRaw(:,i,j),330.46-0,16606,-21.042,101.7686);
end
end

SpectralDataRaw = real(SpectralDataRaw);

RawforNoise = real(RawforNoise);

temp = Noise*randn(numel(SiteA),32,64);
stdn = std(temp(:,1,1));
SNRraw1D = max(SpectralDataRaw,[],'all')/stdn;
SNRrawT1= max(RawforNoise,[],'all')/std(RawforNoise(end,:,end));

if e~=0
   SNRpcaT1 =  max(SpectralDataRaw,[],'all')/std(SpectralDataRaw(end,:,end));
end

% Generate Seperable Relaxation Kernels
MinRateR1 = 0.3;                % 1/s
MaxRateR1 = 10;              % 1/s

MinRateR2 = 70;              % 1/s
MaxRateR2 = 1700;              % 1/s

NumRatesR1 = 100;              % Number of R1 rates
NumRatesR2 = 200;             % Number of R2 rates
R1Rates = logspace(log10(MinRateR1),log10(MaxRateR1),NumRatesR1); % R1 Vector
R2Rates = logspace(log10(MinRateR2),log10(MaxRateR2),NumRatesR2); % R2 Vector
KernelR1 = 1-2*exp(-transpose(t1)*R1Rates);  % R1 Relaxation Kernel
KernelR2 = exp(-transpose(t2)*R2Rates);    % R2 Relaxation Kernel

tic
%% Data Processing
% SVD of Kernels
[Compressed.U_R1, Compressed.S_R1, Compressed.V_R1] = svd(KernelR1);
[Compressed.U_R2, Compressed.S_R2, Compressed.V_R2] = svd(KernelR2);

% Keep 10% of the singular values and compress SVD output matrices
Compressed.CompressionFactor = 0.2; % Percentage of singular values to keep
Compressed.NumberSingularValuesR1 = round(Compressed.CompressionFactor*max(size(diag(Compressed.S_R1))));
Compressed.NumberSingularValuesR2 = round(Compressed.CompressionFactor*max(size(diag(Compressed.S_R2))));

% Compress Sigma
Compressed.S_R1_Truncated = Compressed.S_R1(1:Compressed.NumberSingularValuesR1,1:Compressed.NumberSingularValuesR1);
Compressed.S_R2_Truncated = Compressed.S_R2(1:Compressed.NumberSingularValuesR2,1:Compressed.NumberSingularValuesR2);

% Compress U
Compressed.U_R1_Truncated = Compressed.U_R1(:,1:Compressed.NumberSingularValuesR1);
Compressed.U_R2_Truncated = Compressed.U_R2(:,1:Compressed.NumberSingularValuesR2);

% Compress V
Compressed.V_R1_Truncated = Compressed.V_R1(:,1:Compressed.NumberSingularValuesR1);
Compressed.V_R2_Truncated = Compressed.V_R2(:,1:Compressed.NumberSingularValuesR2);

% Compress Kernels
Compressed.KernelR1_Truncated = Compressed.S_R1_Truncated*transpose(Compressed.V_R1_Truncated);
Compressed.KernelR2_Truncated = Compressed.S_R2_Truncated*transpose(Compressed.V_R2_Truncated);
Kernel = kron(Compressed.KernelR2_Truncated,Compressed.KernelR1_Truncated);

kernel_L1 = get_l(NumRatesR1*NumRatesR2,2);                  
kernel_L2 = get_l(NumRatesR1*NumRatesR2,2);
Kernel_MxyL =[Kernel;alpha*kernel_L1;alpha*kernel_L2];   % Tikhonov regularization with L=2nd discrete deriv. op.

if lambda == 0
    parfor i=1:numel(SiteA)
    % Compress measured data
    DataRaw_Truncated = transpose(Compressed.U_R1_Truncated)*squeeze(SpectralDataRaw(i,:,:))*Compressed.U_R2_Truncated;

    % Compressed Data Stacking 
    DataRawStacked = reshape(DataRaw_Truncated,[size(DataRaw_Truncated,2)*size(DataRaw_Truncated,1),1]);
    DataRawConc = [ DataRawStacked; zeros(size(kernel_L1,1),1); zeros(size(kernel_L1,1),1)];
    %Perform NNLS Fitting of Stacked Data
    %RatesDistribution = lsqnonneg(Kernel,DataRawStacked);
    RatesDistribution = lsqnonneg(Kernel_MxyL,DataRawConc);
    %RatesDistribution = l1_ls_nonneg(Kernel_MxyL,DataRawConc,lambda,1e-3,true);
    % Unstacked the NNLS Output
    RatesDistributionUnStacked(i,:,:) = reshape(RatesDistribution,NumRatesR1,[]);
    disp(i)
    end
else
    parfor i=1:numel(SiteA)
    % Compress measured data
    DataRaw_Truncated = transpose(Compressed.U_R1_Truncated)*squeeze(SpectralDataRaw(i,:,:))*Compressed.U_R2_Truncated;

    % Compressed Data Stacking 
    DataRawStacked = reshape(DataRaw_Truncated,[size(DataRaw_Truncated,2)*size(DataRaw_Truncated,1),1]);
    DataRawConc = [ DataRawStacked; zeros(size(kernel_L1,1),1); zeros(size(kernel_L1,1),1)];
    %Perform NNLS Fitting of Stacked Data
    %RatesDistribution = lsqnonneg(Kernel,DataRawStacked);
    %RatesDistribution = lsqnonneg(Kernel_MxyL,DataRawConc);
    RatesDistribution = l1_ls_nonneg(Kernel_MxyL,DataRawConc,lambda,1e-3,true);
    % Unstacked the NNLS Output
    RatesDistributionUnStacked(i,:,:) = reshape(RatesDistribution,NumRatesR1,[]);
    disp(i)
    end
end
disp('Finished!')
% Axis
sw = 1000;
Freq = [-sw/2:sw/numel(SiteA):sw/2-sw/numel(SiteA)];
Time = toc;

Locate_Spec_Sim_2
