%% Test spike time generator
confi.segment_size  = 3;
confi.segment_shift = 1;
neuronNum_n         = 20;
outputsynNum_k      = 10;
dat2D               = zeros(10,10);
for i=1:100
   dat2D(i) = i; 
end
param.binNumber     = 1000;
param.maxPercSpikes = 0.1;
param.tau           = 0.2;

[times_s, index] = mat2inSpike.poisson(dat2D,param);
subplot(4,1,1)
plot(times_s, index, '.')
[times_s, index] = mat2inSpike.poisson_first(dat2D,param);
subplot(4,1,2)
plot(times_s, index, '.')
[times_s, index] = mat2inSpike.deterministic(dat2D,param);
subplot(4,1,3)
plot(times_s, index, '.')
[times_s, index] = mat2inSpike.expDec(dat2D,param);
subplot(4,1,4)
plot(times_s, index, '.')

%% test convolution
im2D  = imread('cameraman.tif');
[imT] = Dat2Dconvolver.get_dat2Dconv(im2D, 5,'centerON');
figure(1)
subplot(1,2,1);
image(im2D);
subplot(1,2,2);
image(imT);

%% show generated input spikes
im2D  = imread('cameraman.tif');
[imT] = Dat2Dconvolver.get_dat2Dconv(im2D, 5,'centerON');
size_imT = size(imT);
% adapt size to fit segment size and shift
segment_size_data2D  = 4;
segment_shift_data2D = 2;
reduceSize = mod(size_imT(1)-segment_size_data2D,segment_shift_data2D);
imT(end-reduceSize+1:end,:) = [];
imT(:,end-reduceSize+1:end) = [];

param.binNumber     = 1000;
param.maxPercSpikes = 1;
param.tau           = 0.1;

[times_s, index] = mat2inSpike.deterministic(imT, param);

% plot time series
figure;
plot(times_s,index,'.')

%% plot 2D input spikes as movie
fh = figure;
[iR, iC] = size(imT);
axis([1 iR 1 iC], 'manual');
ah = gca;

[matIdxr,matIdxc ]= ind2sub([iR iC], index);
L = length(index);
L_s = 10;
st  = 100; 

Lend = L;

for i=1:st:Lend
    sel = min(1:i,Lend);
    plot(ah, matIdxc(sel),iR-matIdxr(sel)+1,'.b'); hold on;
    if(i>st)
        selNew = i-st:i;
        plot(ah, matIdxc(selNew),iR-matIdxr(selNew)+1,'.r');
    end
    xlim([1 iC]);
    ylim([1 iR]);
    pause(st/Lend * L_s);
end
