%% test convolution
im2D  = imread('cameraman.tif');
[imT] = Dat2Dconvolver.get_dat2Dconv(im2D, 3,'centerON');
figure(1)
subplot(1,2,1);
image(im2D);
subplot(1,2,2);
image(imT);
