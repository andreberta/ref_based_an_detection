close all
clear
clc

train
save 'train.mat'
findThreshold_perSigma
save('thresh.mat','thres_fov_sigma','thres_NL_sigma');