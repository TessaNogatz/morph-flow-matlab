%------------------------------------------------------------------------------
% Usage of morph-flow
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%-----------------------------------------------------------------------------
%
%% Morphological Wavelet Decomposition needs to be available in search path
addpath(genpath('..\morphological-wavelets-matlab'));
addpath(genpath(pwd));

%% Load Sample
load mristack.mat
fixed = im2double(squeeze(mristack));

%% Generate moving by rotating inner part
local = fixed(160:200,100:140,:);

[sx,sy,sz] = deal(0.06,0.02,0.01);
[tx,ty,tz] = deal(0.7,1.2,-0.5);
A = [1 sx 0 tx; sy 1 0 ty; 0 sz 1 tz; 0 0 0 1];
tform = affinetform3d(A);
local = imwarp(local,tform, 'OutputView', imref3d(size(local)));
moving = fixed;
moving(160:200,100:140,:) = local;

fixed = padarray(fixed, [12 12 12], 0, 'both');
moving = deformImg(fixed, 5);%padarray(moving, [12 12 12], 0, 'both');


%% Fix parameters. Note: MorphFlow performs significantly better if 
%% image values are scaled between 0 and 1


[p,q,r] = size(fixed);

para.levels = 12;
para.end_lvl = 0+1;
para.tau = 0.125;
para.lambda = 10;
para.theta = 0.02;
para.warps = 15;
para.epsilon = 0.0001;
para.max_it = 30;
para.method = 'dual';
para.mu =1;
para.gamma = 20.0;
para.outer_iter = 20;
para.min_iter = 5;
para.inner_iter = 5;
para.outer_iter_u = 5;
para.inner_iter_u = 5;
para.median_filtering=1;  %median filtering
para.medianx=5;
para.mediany=5;
para.medianz=5;
para.verbose = 0;


%% Compute morph-flow with given parameters and time sample
tic;

[C0, S0] = QLiftDec3MinMaxMin(fixed,para.levels);
[C1, S1] = QLiftDec3MinMaxMin(moving,para.levels);
    
[u, v, w] = morph_tvl1_of(C0, S0, C1, S1, para);
toc;

%% Show results

figure, imshowpair(fixed(:,:,20),moving(:,:,20))
res0 = abs(fixed-moving);
figure, imagesc(res0(:,:,20), [0 0.99]);

moving_reg = warpImage3D(moving, u, v, w);

res1 = abs(fixed-moving_reg);
figure, imagesc(res1(:,:,20), [0 0.99]);

figure, imshowpair(fixed(:,:,20),moving_reg(:,:,20))
