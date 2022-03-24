% Experimental 3D RAS example for 35Cl Xyla : Isox mixture
% Without regularization this runs in 123 sec. w/ 12 parrallel cores
% Including l2 norm (alpha) will be on the order of minutes
% Including l1 norm (lambda) will be on the order of hours
% ARA + MJJ 2021
clear;
load('3Ddespec.mat') %Load pre-processed 3D spectral data (freq x T1 x T2)
spec = real(despec(198:300,:,:)); %Truncate baseline frequency data
spec = spec /  max(spec, [], 'All');
[a b c] = size(spec);
load('t2Cl.mat');  %load vector of T2 echo times
Tx = 0; %Transmitter in kHz (can adjust later)

%%
t1 = logspace(log10(0.001),log10(3), 32); %Adjust tau vector (VDlist)

alpha = 0; %l2 norm ( = 0.0001 in paper)
lambda = 10; %l1 norm ( = 5 in paper. Long run time!)

% Generate Seperable Relaxation Kernels
MinRateR1 = 0.5;               % 1/s
MaxRateR1 = 1000;              % 1/s

MinRateR2 = 10;              % 1/s
MaxRateR2 = 500;             % 1/s

NumRatesR1 = 300;              % Number of R1 rates
NumRatesR2 = 300;             % Number of R2 rates

R1Rates = logspace(log10(MinRateR1),log10(MaxRateR1),NumRatesR1); % R1 Vector
R2Rates = logspace(log10(MinRateR2),log10(MaxRateR2),NumRatesR2); % R2 Vector

KernelR1 = 1-2*exp(-transpose(t1)*R1Rates); % R1 Relaxation Kernel (can adjust for satrec*)
KernelR2 = exp(-transpose(t2)*R2Rates);    % R2 Relaxation Kernel
tic
%% Data Processing
% SVD of Kernels
[Compressed.U_R1, Compressed.S_R1, Compressed.V_R1] = svd(KernelR1);
[Compressed.U_R2, Compressed.S_R2, Compressed.V_R2] = svd(KernelR2);

%%%Maximum entropy r selection%%%
[m , n] = size(KernelR1);
l = min([m,n]);
s = svd(KernelR1);

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

[notr, Kernelr1 ] = min( delE - (c+d));
%%%Done entropy section - use r for TSVD%%%
%%%Maximum entropy r selection%%%
[m , n] = size(KernelR2);
l = min([m,n]);
s = svd(KernelR2);

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

[notr, Kernelr2 ] = min( delE - (c+d));
%%%Done entropy section - use r for TSVD%%%
Compressed.NumberSingularValuesR1 = Kernelr1;
Compressed.NumberSingularValuesR2 = Kernelr2;

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
Kernel_MxyL =[Kernel;alpha*kernel_L1]; % Tikhonov regularization with L=2nd discrete deriv. op.

if lambda == 0
    parfor i=1:a
    % Compress measured data
    DataRaw_Truncated = transpose(Compressed.U_R1_Truncated)*squeeze(spec(i,:,:))*Compressed.U_R2_Truncated;

    % Compressed Data Stacking 
    DataRawStacked = reshape(DataRaw_Truncated,[size(DataRaw_Truncated,2)*size(DataRaw_Truncated,1),1]);
    DataRawConc = [ DataRawStacked; zeros(size(kernel_L1,1),1)];
    %Perform NNLS Fitting of Stacked Data
    %RatesDistribution = lsqnonneg(Kernel,DataRawStacked);
    RatesDistribution = lsqnonneg(Kernel_MxyL,DataRawConc);
    %RatesDistribution = l1_ls_nonneg(Kernel_MxyL,DataRawConc,lambda,1e-3,true);
    % Unstacked the NNLS Output
    RatesDistributionUnStacked(i,:,:) = reshape(RatesDistribution,NumRatesR1,[]);
    
    disp(i)
    end
else
    parfor i=1:a
    % Compress measured data
    DataRaw_Truncated = transpose(Compressed.U_R1_Truncated)*squeeze(spec(i,:,:))*Compressed.U_R2_Truncated;

    % Compressed Data Stacking 
    DataRawStacked = reshape(DataRaw_Truncated,[size(DataRaw_Truncated,2)*size(DataRaw_Truncated,1),1]);
    DataRawConc = [ DataRawStacked; zeros(size(kernel_L1,1),1)];
    %Perform NNLS Fitting of Stacked Data
    %RatesDistribution = lsqnonneg(Kernel,DataRawStacked);
    %RatesDistribution = lsqnonneg(Kernel_MxyL,DataRawConc);
    RatesDistribution = l1_ls_nonneg(Kernel_MxyL,DataRawConc,lambda,1e-3,true);
    % Unstacked the NNLS Output
    RatesDistributionUnStacked(i,:,:) = reshape(RatesDistribution,NumRatesR1,[]);
    disp(i)
    end
end
disp('All Done!');
Time = toc;
% Axis
sw = 250;
Freq = [-sw/2+Tx:sw/a:sw/2-sw/a+Tx];
Locate_Spec_Xyla_Isox
