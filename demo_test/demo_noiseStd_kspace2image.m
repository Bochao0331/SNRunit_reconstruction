%% This script tests the change of noise Std  during transformation from k-space to image space:
Mat_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0524VOL11Lung\MAT'; % Direcoty where MAT saved
Mat_file = fullfile(Mat_folder,'tse_ES_BL_tra_bh_128_te15_esp15.mat'); % MAT name
load(Mat_file);

%Generate white noise (like we have after prewhitening)
noise_white = complex(randn(size(kdata_all{1}(:,:,1,:))),randn(size(kdata_all{1}(:,:,1,:))));

%Transform from k-space to image: 
noise_test = fft2c(noise_white);
sd1 = std(real(noise_white(:)))
sd2 = std(real(noise_test(:))) % see if sd is still 1?