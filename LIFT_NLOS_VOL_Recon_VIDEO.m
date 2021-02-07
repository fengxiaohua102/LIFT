% ************* This script perform LIFT reconstruciton for NLOS video *******************
% ****** It performs the reconstruction in a volumetric (all frames) manner ********
% ****** the (x,y,t) datacube is not unwarped here ***********
close all;
clear all;
clc;
top_DIR = '.\NLOS\Video\';  % video dataset need refocus
calib_DIR = top_DIR;
filename = 'TwoSquareMovie';
lambda_TV = 4e-1;   % the hyperparameter for the FISTA reconstruction
IMG_INV = true;  % use image relay where image is inverted
% NAVG = 3;
% ******
% 1: Load the calibration data results
% ******
Calib_Res = load([calib_DIR,'Calib_Res.mat']);
Calib_Res = Calib_Res.Calib_Res;
streak_tform = Calib_Res.streak_tform;

% ******
% 2: Load video data 
image_t_video  = load([top_DIR filename '\' filename '.mat']);
image_t = image_t_video.image_t_video;

%% 3: setup solver
u = (1:7)-4; s = 4.5;
sub_img_cnt = round(Calib_Res.sub_img_cnt + u.*s);
options = []; options.INVERT = IMG_INV; options.CROP = true;
options.Deconv = false;  options.Normalize = true;
options.Refocus = true; options.USE_TV = false;
options.sub_img_cnt = round(sub_img_cnt);
% options.sub_img_cnt = round(Calib_Res.cntx_depth(:,2)).'; 

%% nlos parameters
dt = 2e-12;  % 2 ps time resolution
DeltaT = 3e8*dt;
laser_pos = Calib_Res.laser_pos;  % laser pos on the wall, relative to global coo.
grid_pos = Calib_Res.grid_pos; 
    
tic
figure(1);
for K = 1: size(image_t,3)    
    % ******
    % Correct for the streak distortion 
    image_t_frame = imwarp(image_t(:,:,K), streak_tform, 'OutputView', imref2d(size(image_t(:,:,1))));
    
    backG = mean(reshape(image_t_frame(900:1000,200:1000),[],1));
    image_t_frame = image_t_frame-backG;
    image_t_frame(image_t_frame<0) = 0;
    image_t_frame(1:130,:) = 0; 
    image_t_frame = image_t_frame(1:880,:);

    % call the solver
    im_crop = fx_LIFT_ReconVOL(Calib_Res, image_t_frame, lambda_TV, options);
    im_crop = gather(im_crop);
    
    t_form = Calib_Res.H_slant2Plane; t_form = inv(t_form).';
    t_form2 = projective2d(t_form);
    im_crop = warpPerspective(im_crop, t_form2);
    nlos_data = nlos_unwarp(im_crop, Calib_Res, dt);  
    save([top_DIR filename '\NLOS_' filename '_' num2str(K) '.mat'], 'nlos_data', ...
        'DeltaT','grid_pos','laser_pos');
    
    subplot(1,2,1); imagesc(image_t_frame); colormap('hot'); title(['Frame' num2str(K)]); 
    subplot(1,2,2); imagesc(sum(im_crop,3)); colormap('hot'); title('DC image');
    drawnow();
end

toc
function im_data = warpPerspective(im_data, H_slant2Plane)
    ref2d = imref2d(size(im_data(:,:,1)));
    for P = 1:size(im_data,3)
        im_data(:,:,P) = imwarp(im_data(:,:,P), H_slant2Plane, 'OutputView',ref2d);
    end
end