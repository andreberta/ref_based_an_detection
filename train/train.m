%% Prepare workspace
close all
clear
clc

%% Prepare images
%load images
path = 'train_images/';
images = {'fig1.jpg','fig2.jpg','fig3.jpg','fig4.jpg','fig5.jpg',...
    'fig6.png','fig7.jpg','fig8.jpg','fig9.jpg','fig10.gif','fig11.jpg',...
    'fig12.jpg','fig13.jpg','fig14.jpg','fig15.jpg','fig16.jpg'...
    ,'fig17.jpg','fig18.jpg','fig19.jpg','fig20.jpg','fig21.jpg'...
    ,'fig22.jpg','fig23.jpg','fig24.jpg','fig25.jpg','fig26.png'};

images_length = length(images);
images_loaded = cell(1,images_length);

for ii=1:images_length
    im = imread(strcat(path,images{ii}));
    if(ndims(im) == 3)
        im = rgb2gray(im);
    end
    images_loaded{ii} = double(im);
end

%create source images introducing anomaly manually
[rr, cc] = meshgrid(1:150);
A = sqrt((rr-70).^2+(cc-57).^2)<=3;
compare = A;
A = conv2(double(A), ones(5)/5^2, 'same');

B = sqrt((rr-100).^2+(cc-20).^2)<=3;
compare = logical(compare + B);
B = conv2(double(B), ones(6)/5^2, 'same');


images_src = cell(1,images_length);



create_src

%% process images
%prepare data structure for storing data during the computation
sigma = {10,20,30,50,70};
sigma_lenght = length(sigma);

res_NL = cell(sigma_lenght,images_length);
res_fov = cell(sigma_lenght,images_length);
anomaly_NL = cell(sigma_lenght,images_length);
anomaly_fov = cell(sigma_lenght,images_length);

PSNR_NL = cell(sigma_lenght,images_length);
PSNR_fov = cell(sigma_lenght,images_length);
AUC_NL = cell(sigma_lenght,images_length);
AUC_fov = cell(sigma_lenght,images_length);

computation_time_NL = cell(sigma_lenght,images_length);
computation_time_fov = cell(sigma_lenght,images_length);

window_rad_fov = 8;
patch_rad_fov = {3,5,6,8,9};

window_rad_NL = {5,4,4,5,6};
patch_rad_NL = {2,5,5,6,7};


%run for every images in the train set and for every noise level
for jj=1:sigma_lenght
    for ii=1:images_length
        %add noise
        im_ref_noised = images_loaded{ii} + sigma{jj}*randn(size(images_loaded{ii}));
        im_src_noised = images_src{ii} + sigma{jj}*randn(size(images_src{ii}));
        
        %execute the NL-implementation and store computation time
        tstart_NL = tic;
        [res,anomaly] = ...
            NL_impl(im_ref_noised,im_src_noised,window_rad_NL{jj},patch_rad_NL{jj},sigma{jj});
        computation_time_NL{jj,ii} = toc(tstart_NL);
        
        %execute the fov-implementation and store computation time
        tstart_fov = tic;
        [y_hat,patch_hat]=...
            fov_impl(im_ref_noised,im_src_noised,window_rad_fov,patch_rad_fov{jj},sigma{jj});
        computation_time_fov{jj,ii} = toc(tstart_fov);
        
        %store result for NL
        res_NL{jj,ii} = res;          %<-- the reconstructed image
        anomaly_NL{jj,ii} = anomaly;  %<-- the patch wise difference
        
        %store result for fov
        res_fov{jj,ii} = y_hat;          %<-- the reconstructed image
        anomaly_fov{jj,ii} = patch_hat;  %<-- the patch wise difference
        
        %compute and store the reconstruction quality
        PSNR_fov{jj,ii} =10*log10...
            (255^2/mean((images_src{ii}(:)-y_hat(:)).^2));
        PSNR_NL{jj,ii} =10*log10...
            (255^2/mean((images_src{ii}(:)-res(:)).^2));
        
        %compute and store the AUC
        [AUC_NL{jj,ii},~,~] = ROC( anomaly_NL{jj,ii} , compare , 0 );
        [AUC_fov{jj,ii},~,~] = ROC( anomaly_fov{jj,ii} , compare , 0 );
    end
    disp(jj);
end