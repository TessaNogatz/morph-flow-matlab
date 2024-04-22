function [u, v, w] = morph_tvl1_of(C0, S0, C1, S1, para)
%------------------------------------------------------------------------------
% Calculates primal dual optical flow in 3D with morphological coarse to
% fine scheme
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
nx = round(S0(end,4));
ny = round(S0(end,5));
nz = round(S0(end,6));

u = zeros(nx, ny, nz);
v = zeros(nx, ny, nz);
w = zeros(nx, ny, nz);

p1x = zeros(nx, ny, nz); p1y = zeros(nx, ny, nz); p1z = zeros(nx, ny, nz);
p2x = zeros(nx, ny, nz); p2y = zeros(nx, ny, nz); p2z = zeros(nx, ny, nz);
p3x = zeros(nx, ny, nz); p3y = zeros(nx, ny, nz); p3z = zeros(nx, ny, nz);



for l=para.levels:-3:para.end_lvl
    
    % Cartesian voxel setting
    if l >= para.levels
        I0_A000 = retrieveR3D(l, 'a', C0, S0);    
        I1_A000 = retrieveR3D(l, 'a', C1, S1);
        [p, q, r] = size(I0_A000);
        
        u = imresize3(u,  [p, q, r]);
        v = imresize3(v,  [p, q, r]);
        w = imresize3(w,  [p, q, r]);
   
    end
           
    [u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ...
        dual_tvl1_of(I0_A000, I1_A000, ...
        u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,...
        para.tau, para.lambda,para.theta, para.warps, para.epsilon, para.max_it,...
        para.verbose);
    
    tp = class(I0_A000);
    I0_D101 = retrieveR3D(l, 'd', C0, S0);       % "odd slots"
    I1_D101 = retrieveR3D(l, 'd', C1, S1);
    minO = min(min(I0_A000, [], 'all'), min(I0_D101, [], 'all'));
    maxO = max(max(I0_A000, [], 'all'), max(I0_D101, [], 'all'));
    cmin = minO-(maxO-minO);
    cmax = maxO+(maxO-minO);
    
    % reverse update
    sizeA000 = size(I0_A000);
    I0_A000 = I0_A000 - min(zeros(sizeA000,tp),...
          d2Lift000min(I0_D101, size(I0_A000, 1), size(I0_A000, 3), cmax));
    I1_A000 = I1_A000 - min(zeros(sizeA000,tp),...
          d2Lift000min(I1_D101, size(I0_A000, 1), size(I0_A000, 3), cmax));
    
    % reverse predict
    I0_A101 = I0_D101 + ...
           d2Lift101min(I0_A000, size(I0_D101, 1), size(I0_D101, 3), cmax);
    I1_A101 = I1_D101 + ...
           d2Lift101min(I1_A000, size(I0_D101, 1), size(I0_D101, 3), cmax);
    
    % I0_A000 and I0_A101 now are the components of the rectangular voxel
    % setting. 
    
    % Rectangular voxel setting
    % Note that the vector field components are lifted inside the function
    [u_A000, u_A101, v_A000, v_A101, w_A000, w_A101, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ...
        morph_dual_tvl1_of_rect(I0_A000, I0_A101, I1_A000, I1_A101, ...
        u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,...
        para.tau, para.lambda,para.theta, para.warps, para.epsilon, para.max_it, ...
        para.verbose);
    
    [I0_D110, I0_D011] = retrieved1Lift(l-1, 'd', C0, S0);
    [I1_D110, I1_D011] = retrieved1Lift(l-1, 'd', C1, S1);
    
    I0_A000 = I0_A000 - ...
         max(zeros(size(I0_A000),tp), d1Lift000max(I0_D110, I0_D011, size(I0_A000, 2), cmin));
    I0_A101 = I0_A101 - ...
         max(zeros(size(I0_A101),tp), d1Lift101max(I0_D110, I0_D011, size(I0_A101, 2), cmin));
     
    I1_A000 = I1_A000 - ...
         max(zeros(size(I1_A000),tp), d1Lift000max(I1_D110, I1_D011, size(I1_A000, 2), cmin));
    I1_A101 = I1_A101 - ...
         max(zeros(size(I1_A101),tp), d1Lift101max(I1_D110, I1_D011, size(I1_A101, 2), cmin));
     
    I0_A110 = I0_D110 + d1Lift110max(I0_A000, I0_A101, size(I0_D110, 2), cmin);
    I0_A011 = I0_D011 + d1Lift011max(I0_A000, I0_A101, size(I0_D011, 2) ,cmin);
   
    I1_A110 = I1_D110 + d1Lift110max(I1_A000, I1_A101, size(I1_D110, 2), cmin);
    I1_A011 = I1_D011 + d1Lift011max(I1_A000, I1_A101, size(I1_D011, 2) ,cmin);
    
    % Dodecaeder setting. Lifting to dodecaeder and from dodecader back to
    % cartesian for next iteration will be done within the function
    [u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = morph_dual_tvl1_of_dode(...
        I0_A000, I0_A110, I0_A011, I0_A101,...
        I1_A000, I1_A110, I1_A011, I1_A101, ...
        u_A000, u_A101, v_A000, v_A101, w_A000, w_A101,...
        p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,...
        para.tau, para.lambda,para.theta, para.warps, para.epsilon, para.max_it, ....
        para.verbose);
   
     I0_A000 = retrieveImg_hv(C0, S0, l, I0_A000, I0_A110, I0_A011, I0_A101, cmax, 'min');
     I1_A000 = retrieveImg_hv(C1, S1, l, I1_A000, I1_A110, I1_A011, I1_A101, cmax, 'min');
     if l-2 == para.end_lvl
        % interpolate to get the initial value of the finner pyramid level
       [u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ...
        dual_tvl1_of(I0_A000, I1_A000, ...
        u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,...
        para.tau, para.lambda,para.theta, para.warps, para.epsilon, para.max_it,...
        para.verbose);
     end
end