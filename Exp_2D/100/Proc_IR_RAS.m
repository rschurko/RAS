clear;
data=bread('ser');
%Parameters to adjust:
dw = 2;
np = 300/dw;              % Number of points in 1 echo
np_dead = (50 + 10 +2)/dw; % Number of dead points between echoes
nechoes = 87;         % Number of spin echoes
zf = 512;             %Number of points to be processed
nt = 32;

tau = logspace(log10(0.1),log10(12),32); %Adjust tau vector (VDList)

data = reshape(data,numel(data)/nt,nt);

for i = 1:nt
    fid = data(:,i);
    %Chop of DF
    fid = transpose(fid);
    DF_size = 123;            % DECIMATION SIZE - THIS MAY NEED TO BE CHANGED
    fid = fid(DF_size:end);

    %Eliminate dead time spacing
    for n = 1:nechoes
        fid(np*n:(np*n + (np_dead-1))) = [];
    end
    %Eliminate trailing dead time points
    fid((np*nechoes + 1):end) = [];
    
        %sd is like inverse of lb
sd = 100;
for z = 1:1:np
    n = z - np/2;
    gauss(z) = ((1/(2*3.14*sd))*exp(-((n-1)^2)/(2*sd^2)));
end

    %Create CPMG matrix of fid
    CPMG = reshape(fid,np,nechoes);
    for z=1:1:nechoes
      CPMG(:,z) = CPMG(:,z).*gauss';
    end
    CPMG = [CPMG; zeros((zf-np),nechoes)];
for j = 1:nechoes
    sp(:,j) = fftshift(fft(CPMG(:,j)));
    %sp2 = phi(sp,55,29000,-9700);
end
    raw(:,i,:) = sp;
    sp2 = phi_v3(sp,158.0-180, 18313.7, -2098.2, -17);
    Z = denoise(sp,3);
    sp3 = phi_v3(Z,158.0-180, 18313.7, -2098.2, -17);
    spec(:,i,:) = sp2;
    %area(i) = trapz(spec(441:592,i));
    despec(:,i,:) = sp3;
end

%%Std noise measure
nspec = spec / max(spec, [], 'All');
ndespec = despec / max(despec, [], 'All');

rawstd = std(real(nspec(1:100,1,2)));
destd = std(real(ndespec(1:100,1,2)));
%%end

for i = 1:nechoes
    rawt1 = raw(:,:,i);
    Z1 = denoise(rawt1,4);
    sp4 = phi_v3(Z1,158.0-180, 18313.7, -2098.2, -17);    
    despect1(:,:,i) = sp4;
end

subplot(2,2,1)
mesh(real(spec(:,:,2)));
subplot(2,2,2)
mesh(real(despect1(:,:,2)));
subplot(2,2,3)
mesh(real(squeeze(spec(:,end,:))))
subplot(2,2,4)
mesh(real(squeeze(despect1(:,end,:))))

rawspecR1 = real(sum(spec,3));
rawspecR2 = squeeze(real(sum(abs(spec),2)));

despecR1 = real(sum(despect1,3));
despecR2 = squeeze(real(sum(abs(despect1),2)));

for k=1:150
SNRrawT1(k) = max(real(nspec),[],'all')/std(real(nspec(315+k,:,1)));
SNRpcaT1(k) =  max(real(ndespec),[],'all')/std(real(ndespec(315+k,:,1)));
SNRrawT2(k) = max(real(nspec),[],'all')/std(real(nspec(315+k,end,:)));
SNRpcaT2(k) = max(real(ndespec),[],'all')/std(real(ndespec(315+k,end,:)));
end

SNRrawT1 = mean(SNRrawT1);
SNRpcaT1 = mean(SNRpcaT1);
SNRrawT2 = mean(SNRrawT2);
SNRpcaT2 = mean(SNRpcaT2);

disp('T1 SNR epsilon ='); disp(SNRpcaT1/SNRrawT1);
disp('T2 SNR epsilon ='); disp(SNRpcaT2/SNRrawT2);