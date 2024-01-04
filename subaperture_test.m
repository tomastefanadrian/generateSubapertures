clearvars
clc
close all

load('HV.mat')
a=0.75;

%38-7+38-7+38
%30-7+30-7+30
winSize=30;
overlap=7;
N=100;

imHH=reshape(HHHV(1,:,:),100,100);
imHV=reshape(HHHV(2,:,:),100,100);

figure
imagesc(abs(imHV)),colormap('Gray')

HH=fft(imHH);
HV=fft(imHV);

w=coefHamming(N,a);

nHH=zeros(size(HH));
nHV=zeros(size(HV));
invW=1./w;
for i=1:100
    nHH(:,i)=invW.*fftshift(HH(:,i));
    nHV(:,i)=invW.*fftshift(HV(:,i));
end

energyProfile=abs(diff(sum(abs(nHH),2)));
valMaxEnergyProfile=max(energyProfile);
[~,locs]=findpeaks(energyProfile,'MinPeakHeight',valMaxEnergyProfile*0.8);
zeroPadding=zeros(locs(1),100);
usefulSpectrumIndex=[locs(1),locs(2)];
usefulSpectrumSize=locs(2)-locs(1)+1;

usefulHH=nHH(usefulSpectrumIndex(1):usefulSpectrumIndex(2),:);
usefulHV=nHV(usefulSpectrumIndex(1):usefulSpectrumIndex(2),:);

newWin=coefHamming(winSize,a);
noWin=ones(winSize,1);

usefulWinHH=usefulHH(1:winSize,:);
for i=1:N
    usefulWinHH(:,i)=newWin.*usefulWinHH(:,i);
end
s1HH=[usefulWinHH;zeros(usefulSpectrumSize-winSize,N)];

usefulWinHH=usefulHH(winSize-overlap+1:winSize-overlap+winSize,:);
for i=1:N
    usefulWinHH(:,i)=newWin.*usefulWinHH(:,i);
end
s2HH=[zeros(winSize-overlap,N);usefulWinHH;zeros(usefulSpectrumSize-2*winSize+overlap,N)];

usefulWinHH=usefulHH(2*(winSize-overlap)+1:2*(winSize-overlap)+winSize,:);
for i=1:N
    usefulWinHH(:,i)=newWin.*usefulWinHH(:,i);
end
s3HH=[zeros(2*(winSize-overlap),N);usefulWinHH];


usefulWinHV=usefulHV(1:winSize,:);
for i=1:N
    usefulWinHV(:,i)=newWin.*usefulWinHV(:,i);
end
s1HV=[usefulWinHV;zeros(usefulSpectrumSize-winSize,N)];

usefulWinHV=usefulHV(winSize-overlap+1:winSize-overlap+winSize,:);
for i=1:N
    usefulWinHV(:,i)=newWin.*usefulWinHV(:,i);
end
s2HV=[zeros(winSize-overlap,N);usefulWinHV;zeros(usefulSpectrumSize-2*winSize+overlap,N)];

usefulWinHV=usefulHV(2*(winSize-overlap)+1:2*(winSize-overlap)+winSize,:);
for i=1:N
    usefulWinHV(:,i)=newWin.*usefulWinHV(:,i);
end
s3HV=[zeros(2*(winSize-overlap),N);usefulWinHV];


% s1HH=[usefulHH(1:winSize,:);zeros(usefulSpectrumSize-winSize,N)];
% s2HH=[zeros(winSize-overlap,N);usefulHH(winSize-overlap+1:winSize-overlap+winSize,:);zeros(usefulSpectrumSize-2*winSize+overlap,N)];
% s3HH=[zeros(2*(winSize-overlap),N);usefulHH(2*(winSize-overlap)+1:2*(winSize-overlap)+winSize,:)];
% 
% s1HV=[usefulHV(1:winSize,:);zeros(usefulSpectrumSize-winSize,N)];
% s2HV=[zeros(winSize-overlap,N);usefulHV(winSize-overlap+1:winSize-overlap+winSize,:);zeros(usefulSpectrumSize-2*winSize+overlap,N)];
% s3HV=[zeros(2*(winSize-overlap),N);usefulHV(2*(winSize-overlap)+1:2*(winSize-overlap)+winSize,:)];

s1HH=[zeroPadding;s1HH;zeroPadding];
s2HH=[zeroPadding;s2HH;zeroPadding];
s3HH=[zeroPadding;s3HH;zeroPadding];

s1HV=[zeroPadding;s1HV;zeroPadding];
s2HV=[zeroPadding;s2HV;zeroPadding];
s3HV=[zeroPadding;s3HV;zeroPadding];

% tmp=[newWin;zeros(N-winSize,1)];
% for i=1:100
%     s1HH=tmp.*s1HH;
% end
% 
% tmp=[zeros(winSize-overlap,1);newWin;zeros(winSize-overlap,1)];
% for i=1:100
%     s2HH=tmp.*s2HH;
% end
% 
% tmp=[zeros(winSize-overlap,1);newWin;zeros(winSize-overlap,1)];
% for i=1:100
%     s3HH=tmp.*s3HH;
% end

imnHH=ifft(nHH);
imnHV=ifft(nHV);

ims1HH=ifft(ifftshift(s1HH,1));
ims2HH=ifft(ifftshift(s2HH,1));
ims3HH=ifft(ifftshift(s3HH,1));

ims1HV=ifft(ifftshift(s1HV,1));
ims2HV=ifft(ifftshift(s2HV,1));
ims3HV=ifft(ifftshift(s3HV,1));


% figure
% imagesc(fftshift(abs(HH))),colormap('Gray')
% 
% figure
% subplot(1,3,1),imagesc(abs(ifftshift(s1HH,1))),colormap('Gray')
% subplot(1,3,2),imagesc(abs(ifftshift(s2HH,1))),colormap('Gray')
% subplot(1,3,3),imagesc(abs(ifftshift(s3HH,1))),colormap('Gray')

figure
subplot(2,2,1),imagesc(abs(ims1HH)),colormap('Gray'), title('S1')
subplot(2,2,2),imagesc(abs(ims2HH)),colormap('Gray'), title('S2')
subplot(2,2,3),imagesc(abs(ims3HH)),colormap('Gray'), title('S3')
subplot(2,2,4),imagesc(abs(imHH)),colormap('Gray'), title('Original')

figure
subplot(2,2,1),imagesc(abs(ims1HV)),colormap('Gray'), title('S1')
subplot(2,2,2),imagesc(abs(ims2HV)),colormap('Gray'), title('S2')
subplot(2,2,3),imagesc(abs(ims3HV)),colormap('Gray'), title('S3')
subplot(2,2,4),imagesc(abs(imHV)),colormap('Gray'), title('Original')
