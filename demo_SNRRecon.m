%   This script does reading, and image reconstruction in SNR units for 2D TSE lung data

%   1. You can read either TWIX or MAT(generated based on TWIX) fiels.

%	2. For ISMRMD data, we need raw data, noise data, info

%   3. MAT file contains cells (each cells contains one echo-shift data; same echo-shift order as TWIX):
%       (1) k-space date: k_data_all (Nx x Ny x Ncoil x Nz)
%       (2) coil-combined images (magnitude,sum-of-square combination): img_List_all (Nx x Ny x Nslice)
%       (3) coil images: c_img_all (Nx x Ny x Ncoil x Nslice)Nslice
%       (4) info: information obtained from TWIX (TE,TR,TurboFactor...)

%   Example Raw data locate in: /mnt/sdata_new/Bochao/0528AppleESSE/RAW/ on MREL
%   server;
%   Example MAT files locate in: /mnt/sdata_new/Bochao/0528AppleESSE//MAT

%	Author: Bochao Li
%	Email: bochaoli@usc.edu

clear
%% add path
addpath(genpath('.\functions'));
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
    data_parent_directory = '';
end
addpath(genpath(ismrmrd_directory));
%% ---- Reading parameter settings ----
zero_shift_ind = 1; % zero-shift data index
data_format = 'mat'; % Data reading option 'ismrmd'(Raw) or 'mat';
Mat_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0524VOL11Lung\MAT'; % Direcoty where MAT saved
Mat_file = fullfile(Mat_folder,'tse_ES_BL_tra_bh_128_te11_esp11.mat'); % MAT name
slice_index = 1; % sincle index
%% ---- Loading data ----
if strcmp(data_format, 'ismrmd')
    h5_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0524VOL11Lung\RAW\h5';
    noise_folder = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0524VOL11Lung\RAW\noise';
    h5_fileList = fullfile(h5_folder, 'meas_MID00032_FID05399_tse_tra_bh_128_te11_esp11.h5'); % Change to whatever pattern you need.
    noise_fileList = fullfile(noise_folder, 'noise_meas_MID00032_FID05399_tse_tra_bh_128_te11_esp11.h5'); % Change to whatever pattern you need.
    h5_Files = dir(h5_fileList);
    noise_Files = dir(noise_fileList);
    nFiles = length(h5_Files);
    
    % ----ISMRMD .h5 reading flags----
    Read_flags.h5_fileList       = h5_fileList;
    Read_flags.noise_fileList    = noise_fileList;
    Read_flags.RemoveOS          = true; % remove oversampling
    Read_flags.IgnoreSeg         = true; % concatanate segmemtns
    Read_flags.DoAverage         = true; % Do averages (if 'average' was set during data acquistion)
    Read_flags.CropPhaseEncoding = true;
    Read_flags.Squeeze           = true;
    Read_flags.os                = 2; % oversampling rate (Siemens Default value, don't change it)
    Read_flags.noisePreWhitening = true;
    
    % ----Loading Raw Data----
    [c_img_all, kdata_all, noise_all, info_all] = readTSE_ismrmd(Read_flags); %loading data from TWIX.
    if ~isfile(Mat_file)
        save (Mat_file, 'c_img_all', 'kdata_all','noise_all','info_all');
    end
    
elseif strcmp(data_format, 'mat')
    load(Mat_file);
end

%% ---- Get data info ----
[Nx, Ny, Nslice, Ncoil ] = size(kdata_all{1})
TE = info_all{1}.TE; % msec
TR = info_all{1}.TR; % msec
echoshifts = [0]; % echo shift value
echo_time_all = TE + echoshifts; % Echo time of the echo-shift

%% ---- Extract data list a specifc echo-shift data for seperate fitting&analysis ----
clear img_List c_img kdata c_img_List
order = 'all'; % p: postive shift; n: negative shift; all: all echo-shift data
if strcmp(order, 'p')
    index = zero_shift_ind : length(kdata_all); % For positive echo-shift
elseif strcmp(order, 'n')
    index = 1 : zero_shift_ind;
elseif strcmp(order, 'all')
    index = 1 : length(kdata_all);
end

for k = 1:length(index)
    idx_k = index(k);
    c_img{k} = c_img_all{idx_k}(:,:,slice_index,:);
    kdata{k} = kdata_all{idx_k}(:,:,slice_index,:);
    echotime(k) = echo_time_all(idx_k);
    noise{k} = noise_all{idx_k};
    info{k} = info_all{idx_k};
end
%% Apply the noise prewhitening matrix on k-space before Recon
% ---- Calculate noise covariance ----
for k = 1:length(index)
    [Psi, inv_L] = calculate_noise_covariance(noise{k});
    kdata_prew{k} = prewhitening(kdata{k},inv_L);
end


%% Coil senstivity estimation and
csm_option.method = 'walsh'; % coil estimation method: walsh, sos
csm_option.cal_shape = [16 16]; %calibration region
csm_option.kdata = kdata_prew{1};
[csm, cal_im]  = coil_estimation(csm_option);

%% Root-sum-of-squares (RSS) reconstruction
p = fft2c(kdata_prew{1}); % Ncoil x Nsample (Nx*Ny*Nslice)
snr_rss = sqrt(2)*sqrt(sum(abs(p).^2,4));
snr_walsh=  sum( conj(reshape(csm,[Nx, Ny, size(kdata_prew,3), Ncoil])) .* fft2c(kdata_prew{1}), 4) ./ sqrt(sum(abs(csm).^2,4)); % Coil combined images

figure,imagesc(abs(rot90(snr_walsh,-1)),[0 20]);
colormap(gray)
axis off image
colorbar

figure,imagesc((rot90(snr_rss,-1)),[0 20]);
colormap(gray)
axis off image
colorbar