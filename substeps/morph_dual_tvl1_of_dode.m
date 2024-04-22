function [u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ...
    morph_dual_tvl1_of_dode(...
    I0_A000, I0_A110, I0_A011, I0_A101,...
    I1_A000, I1_A110, I1_A011, I1_A101, ... ...
    u_A000, u_A101, v_A000, v_A101, w_A000, w_A101,...
    p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, ...
    tau, lambda, theta, warps, epsilon, max_it, verbose)

%------------------------------------------------------------------------------
% Calculates primal dual optical flow in 3D on dodecahedral voxel grids
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
[nx, ny, nz] = size(I0_A000);
N = nx*4*ny*4*nz*4;
l_t = lambda * theta;

minO = min(min(u_A000, [], 'all'), 0);
maxO = max(max(u_A000, [], 'all'), 0);
cmin = minO-(maxO-minO);
cmax = maxO+(maxO-minO);

xD110 = size(I0_A110, 2);

[u_A000, u_A110, u_A011, u_A101] = ...
    lift_d1_nodetail(u_A000, u_A101, xD110, cmin, 'max');
[v_A000, v_A110, v_A011, v_A101] = ...
    lift_d1_nodetail(v_A000, v_A101, xD110, cmin, 'max');
[w_A000, w_A110, w_A011, w_A101] = ...
    lift_d1_nodetail(w_A000, w_A101, xD110, cmin, 'max');


[p1x_A000, p1x_A110, p1x_A011, p1x_A101] = ...
    lift_d1_nodetail(p1x{1}, p1x{2}, xD110, cmin, 'max');
[p1y_A000, p1y_A110, p1y_A011, p1y_A101] = ...
    lift_d1_nodetail(p1y{1}, p1y{2}, xD110, cmin, 'max');
[p1z_A000, p1z_A110, p1z_A011, p1z_A101] = ...
    lift_d1_nodetail(p1z{1}, p1z{2}, xD110, cmin, 'max');

[p2x_A000, p2x_A110, p2x_A011, p2x_A101] = ...
    lift_d1_nodetail(p2x{1}, p2x{2}, xD110, cmin, 'max');
[p2y_A000, p2y_A110, p2y_A011, p2y_A101] = ...
    lift_d1_nodetail(p2y{1}, p2y{2}, xD110, cmin, 'max');
[p2z_A000, p2z_A110, p2z_A011, p2z_A101] = ...
    lift_d1_nodetail(p2z{1}, p2z{2}, xD110, cmin, 'max');

[p3x_A000, p3x_A110, p3x_A011, p3x_A101] = ...
    lift_d1_nodetail(p3x{1}, p3x{2}, xD110, cmin, 'max');
[p3y_A000, p3y_A110, p3y_A011, p3y_A101] = ...
    lift_d1_nodetail(p3y{1}, p3y{2}, xD110, cmin, 'max');
[p3z_A000, p3z_A110, p3z_A011, p3z_A101] = ...
    lift_d1_nodetail(p3z{1}, p3z{2}, xD110, cmin, 'max');
   

for i = 1:warps
    
    I1w_A000 = warpImage3D(I1_A000, u_A000, v_A000, w_A000);
    I1w_A110 = warpImage3D(I1_A110, u_A110, v_A110, w_A110);
    I1w_A011 = warpImage3D(I1_A011, u_A011, v_A011, w_A011);
    I1w_A101 = warpImage3D(I1_A101, u_A101, v_A101, w_A101);
    % compute the derivatives, smoothing considering the directions
    [Ix0_A000, Iy0_A000, Iz0_A000] = ...
    gradient000hv(I0_A000, I0_A110, I0_A011, I0_A101, 0);
    [Ix0_A110, Iy0_A110, Iz0_A110] = ...
    gradient110hv(I0_A110, I0_A000, I0_A011, I0_A101, 0);
    [Ix0_A011, Iy0_A011, Iz0_A011] = ...
    gradient011hv(I0_A011, I0_A000, I0_A110, I0_A101, 0);
    [Ix0_A101, Iy0_A101, Iz0_A101] = ...
    gradient101hv(I0_A101, I0_A000, I0_A110, I0_A011, 0);
    
    [Ix1_A000, Iy1_A000, Iz1_A000] = ...
    gradient000hv(I1_A000, I1_A110, I1_A011, I1_A101, 0);
    [Ix1_A110, Iy1_A110, Iz1_A110] = ...
    gradient110hv(I1_A110, I1_A000, I1_A011, I1_A101, 0);
    [Ix1_A011, Iy1_A011, Iz1_A011] = ...
    gradient011hv(I1_A011, I1_A000, I1_A110, I1_A101, 0);
    [Ix1_A101, Iy1_A101, Iz1_A101] = ...
    gradient101hv(I1_A101, I1_A000, I1_A110, I1_A011, 0);
   
    
    %warp the derivative according to the current displacement
    I1wx_A000 = warpImage3D(Ix1_A000, u_A000, v_A000, w_A000); 
    I1wx_A110 = warpImage3D(Ix1_A110, u_A110, v_A110, w_A110);
    I1wx_A011 = warpImage3D(Ix1_A011, u_A011, v_A011, w_A011); 
    I1wx_A101 = warpImage3D(Ix1_A101, u_A101, v_A101, w_A101);
    
    I1wy_A000 = warpImage3D(Iy1_A000, u_A000, v_A000, w_A000);  
    I1wy_A110 = warpImage3D(Iy1_A110, u_A110, v_A110, w_A110);
    I1wy_A011 = warpImage3D(Iy1_A011, u_A011, v_A011, w_A011); 
    I1wy_A101 = warpImage3D(Iy1_A101, u_A101, v_A101, w_A101);
    
    I1wz_A000 = warpImage3D(Iz1_A000, u_A000, v_A000, w_A000); 
    I1wz_A110 = warpImage3D(Iz1_A110, u_A110, v_A110, w_A110);
    I1wz_A011 = warpImage3D(Iz1_A011, u_A011, v_A011, w_A011);  
    I1wz_A101 = warpImage3D(Iz1_A101, u_A101, v_A101, w_A101);
    
    %estimate derivative direction by mixing i0 and i1 derivatives
    I1wx_A000 = 0.5 * (I1wx_A000 + Ix0_A000); 
    I1wx_A011 = 0.5 * (I1wx_A011 + Ix0_A011); 
    I1wx_A110 = 0.5 * (I1wx_A110 + Ix0_A110); 
    I1wx_A101 = 0.5 * (I1wx_A101 + Ix0_A101); 
    
    I1wy_A000 = 0.5 * (I1wy_A000 + Iy0_A000); 
    I1wy_A011 = 0.5 * (I1wy_A011 + Iy0_A011); 
    I1wy_A110 = 0.5 * (I1wy_A110 + Iy0_A110); 
    I1wy_A101 = 0.5 * (I1wy_A101 + Iy0_A101); 
    
    I1wz_A000 = 0.5 * (I1wz_A000 + Iz0_A000); 
    I1wz_A011 = 0.5 * (I1wz_A011 + Iz0_A011); 
    I1wz_A110 = 0.5 * (I1wz_A110 + Iz0_A110); 
    I1wz_A101 = 0.5 * (I1wz_A101 + Iz0_A101); 
    
    Ix2_A000 = I1wx_A000.*I1wx_A000;
    Ix2_A011 = I1wx_A011.*I1wx_A011;
    Ix2_A110 = I1wx_A110.*I1wx_A110;
    Ix2_A101 = I1wx_A101.*I1wx_A101;
    
    Iy2_A000 = I1wy_A000.*I1wy_A000;
    Iy2_A011 = I1wy_A011.*I1wy_A011;
    Iy2_A110 = I1wy_A110.*I1wy_A110;
    Iy2_A101 = I1wy_A101.*I1wy_A101;
    
    Iz2_A000 = I1wz_A000.*I1wz_A000;
    Iz2_A011 = I1wz_A011.*I1wz_A011;
    Iz2_A110 = I1wz_A110.*I1wz_A110;
    Iz2_A101 = I1wz_A101.*I1wz_A101;
    
    grad_A000 = Ix2_A000 + Iy2_A000 + Iz2_A000;
    grad_A110 = Ix2_A110 + Iy2_A110 + Iz2_A110;
    grad_A011 = Ix2_A011 + Iy2_A011 + Iz2_A011;
    grad_A101 = Ix2_A101 + Iy2_A101 + Iz2_A101;
    
    rho_c_A000 = I1w_A000 - I0_A000 - ...
        I1wx_A000.*u_A000 - I1wy_A000.*v_A000 - I1wz_A000.*w_A000;
    rho_c_A110 = I1w_A110 - I0_A110 - ...
        I1wx_A110.*u_A110 - I1wy_A110.*v_A110 - I1wz_A110.*w_A110;
    rho_c_A011 = I1w_A011 - I0_A011 - ...
        I1wx_A011.*u_A011 - I1wy_A011.*v_A011 - I1wz_A011.*w_A011;
    rho_c_A101 = I1w_A101 - I0_A101 - ...
        I1wx_A101.*u_A101 - I1wy_A101.*v_A101 - I1wz_A101.*w_A101;
    
    n = 0; error = inf;

    mf = 3;

    u_A000 = medfilt3(u_A000,[mf,mf,mf]);
    v_A000 = medfilt3(v_A000,[mf,mf,mf]);
    w_A000 = medfilt3(w_A000,[mf,mf,mf]);
    
    u_A110 = medfilt3(u_A110,[mf,mf,mf]);
    v_A110 = medfilt3(v_A110,[mf,mf,mf]);
    w_A110 = medfilt3(w_A110,[mf,mf,mf]);
    
    u_A011 = medfilt3(u_A011,[mf,mf,mf]);
    v_A011 = medfilt3(v_A011,[mf,mf,mf]);
    w_A011 = medfilt3(w_A011,[mf,mf,mf]);
    
    u_A101 = medfilt3(u_A101,[mf,mf,mf]);
    v_A101 = medfilt3(v_A101,[mf,mf,mf]);
    w_A101 = medfilt3(w_A101,[mf,mf,mf]);

    while error > epsilon && n < max_it
        n = n + 1;
        rho_A000 = rho_c_A000 + ...
            I1wx_A000 .* u_A000 + I1wy_A000 .* v_A000  + I1wz_A000 .* w_A000;
        rho_msk_1 = rho_A000 < -l_t * grad_A000;
        d1_A000 = l_t * rho_msk_1 .* I1wx_A000;
        d2_A000 = l_t * rho_msk_1 .* I1wy_A000;
        d3_A000 = l_t * rho_msk_1 .* I1wz_A000;
        
        rho_msk_2 = rho_A000 > l_t * grad_A000;
        d1_A000 = d1_A000 - l_t * rho_msk_2 .* I1wx_A000;
        d2_A000 = d2_A000 - l_t * rho_msk_2 .* I1wy_A000;
        d3_A000 = d3_A000 - l_t * rho_msk_2 .* I1wz_A000;
        
        rho_msk_else = ones(size(rho_msk_1)) - rho_msk_1 - rho_msk_2;
        if min(rho_msk_else, [], 'all') >= 0
            d1_A000 = d1_A000 - ...
                    rho_msk_else .* rho_A000 ./ grad_A000 .* I1wx_A000;
            d2_A000 = d2_A000 - ...
                    rho_msk_else .* rho_A000 ./ grad_A000 .* I1wy_A000;
            d3_A000 = d3_A000 - ...
                     rho_msk_else .* rho_A000 ./ grad_A000 .* I1wz_A000;
            
            d1_A000(isnan(d1_A000)) = 0;
            d2_A000(isnan(d2_A000)) = 0;
            d3_A000(isnan(d3_A000)) = 0;
            
        else
            sprintf('something went wrong');
        end
        
        rho_A110 = rho_c_A110 + ...
            I1wx_A110 .* u_A110 + I1wy_A110 .* v_A110  + I1wz_A110 .* w_A110;
        rho_msk_1 = rho_A110 < -l_t * grad_A110;
        d1_A110 = l_t * rho_msk_1 .* I1wx_A110;
        d2_A110 = l_t * rho_msk_1 .* I1wy_A110;
        d3_A110 = l_t * rho_msk_1 .* I1wz_A110;
        
        rho_msk_2 = rho_A110 > l_t * grad_A110;
        d1_A110 = d1_A110 - l_t * rho_msk_2 .* I1wx_A110;
        d2_A110 = d2_A110 - l_t * rho_msk_2 .* I1wy_A110;
        d3_A110 = d3_A110 - l_t * rho_msk_2 .* I1wz_A110;
        
        rho_msk_else = ones(size(rho_msk_1)) - rho_msk_1 - rho_msk_2;
        if min(rho_msk_else, [], 'all') >= 0
            d1_A110 = d1_A110 - ...
                rho_msk_else .* rho_A110 ./ grad_A110 .* I1wx_A110;
            d2_A110 = d2_A110 - ...
                rho_msk_else .* rho_A110 ./ grad_A110 .* I1wy_A110;
            d3_A110 = d3_A110 - ...
                rho_msk_else .* rho_A110 ./ grad_A110 .* I1wz_A110;
            
            d1_A110(isnan(d1_A110)) = 0;
            d2_A110(isnan(d2_A110)) = 0;
            d3_A110(isnan(d3_A110)) = 0;
            
        else
            sprintf('something went wrong');
        end
        
        rho_A011 = rho_c_A011 + ...
            I1wx_A011 .* u_A011 + I1wy_A011 .* v_A011  + I1wz_A011 .* w_A011;
        rho_msk_1 = rho_A011 < -l_t * grad_A011;
        d1_A011 = l_t * rho_msk_1 .* I1wx_A011;
        d2_A011 = l_t * rho_msk_1 .* I1wy_A011;
        d3_A011 = l_t * rho_msk_1 .* I1wz_A011;
        
        rho_msk_2 = rho_A011 > l_t * grad_A011;
        d1_A011 = d1_A011 - l_t * rho_msk_2 .* I1wx_A011;
        d2_A011 = d2_A011 - l_t * rho_msk_2 .* I1wy_A011;
        d3_A011 = d3_A011 - l_t * rho_msk_2 .* I1wz_A011;
        
        rho_msk_else = ones(size(rho_msk_1)) - rho_msk_1 - rho_msk_2;
        if min(rho_msk_else, [], 'all') >= 0
            d1_A011 = d1_A011 - ...
                rho_msk_else .* rho_A011 ./ grad_A011 .* I1wx_A011;
            d2_A011 = d2_A011 - ...
                rho_msk_else .* rho_A011 ./ grad_A011 .* I1wy_A011;
            d3_A011 = d3_A011 - ...
                rho_msk_else .* rho_A011 ./ grad_A011 .* I1wz_A011;
            
            d1_A011(isnan(d1_A011)) = 0;
            d2_A011(isnan(d2_A011)) = 0;
            d3_A011(isnan(d3_A011)) = 0;
            
        else
            sprintf('something went wrong');
        end
        
        rho_A101 = rho_c_A101 + ...
            I1wx_A101 .* u_A101 + I1wy_A101 .* v_A101  + I1wz_A101 .* w_A101;
        rho_msk_1 = rho_A101 < -l_t * grad_A101;
        d1_A101 = l_t * rho_msk_1 .* I1wx_A101;
        d2_A101 = l_t * rho_msk_1 .* I1wy_A101;
        d3_A101 = l_t * rho_msk_1 .* I1wz_A101;
        
        rho_msk_2 = rho_A101 > l_t * grad_A101;
        d1_A101 = d1_A101 - l_t * rho_msk_2 .* I1wx_A101;
        d2_A101 = d2_A101 - l_t * rho_msk_2 .* I1wy_A101;
        d3_A101 = d3_A101 - l_t * rho_msk_2 .* I1wz_A101;
        
        rho_msk_else = ones(size(rho_msk_1)) - rho_msk_1 - rho_msk_2;
        if min(rho_msk_else, [], 'all') >= 0
            d1_A101 = d1_A101 - ...
                rho_msk_else .* rho_A101 ./ grad_A101 .* I1wx_A101;
            d2_A101 = d2_A101 - ...
                rho_msk_else .* rho_A101 ./ grad_A101 .* I1wy_A101;
            d3_A101 = d3_A101 - ...
                rho_msk_else .* rho_A101 ./ grad_A101 .* I1wz_A101;
            
            d1_A101(isnan(d1_A101)) = 0;
            d2_A101(isnan(d2_A101)) = 0;
            d3_A101(isnan(d3_A101)) = 0;
            
        else
            sprintf('something went wrong');
        end
         
        u_tilde_A000 = u_A000 + d1_A000;
        v_tilde_A000 = v_A000 + d2_A000;
        w_tilde_A000 = w_A000 + d3_A000;
        
        u_tilde_A110 = u_A110 + d1_A110;
        v_tilde_A110 = v_A110 + d2_A110;
        w_tilde_A110 = w_A110 + d3_A110;
        
        u_tilde_A011 = u_A011 + d1_A011;
        v_tilde_A011 = v_A011 + d2_A011;
        w_tilde_A011 = w_A011 + d3_A011;
        
        u_tilde_A101 = u_A101 + d1_A101;
        v_tilde_A101 = v_A101 + d2_A101;
        w_tilde_A101 = w_A101 + d3_A101;

        [div_p1_A000, div_p1_A110, div_p1_A011, div_p1_A101] = div_p3D_hv(...
            p1x_A000, p1x_A110, p1x_A011, p1x_A101, ...
            p1y_A000, p1y_A110, p1y_A011, p1y_A101, ...
            p1z_A000, p1z_A110, p1z_A011, p1z_A101);
        [div_p2_A000, div_p2_A110, div_p2_A011, div_p2_A101] = div_p3D_hv(...
            p2x_A000, p2x_A110, p2x_A011, p2x_A101, ...
            p2y_A000, p2y_A110, p2y_A011, p2y_A101, ...
            p2z_A000, p2z_A110, p2z_A011, p2z_A101);
        [div_p3_A000, div_p3_A110, div_p3_A011, div_p3_A101] = div_p3D_hv(...
            p3x_A000, p3x_A110, p3x_A011, p3x_A101, ...
            p3y_A000, p3y_A110, p3y_A011, p3y_A101, ...
            p3z_A000, p3z_A110, p3z_A011, p3z_A101);
        
        error = 0.0;
        
        uk_A000 = u_A000;
        vk_A000 = v_A000;
        wk_A000 = w_A000;
        
        uk_A110 = u_A110;
        vk_A110 = v_A110;
        wk_A110 = w_A110;
        
        uk_A011 = u_A011;
        vk_A011 = v_A011;
        wk_A011 = w_A011;
        
        uk_A101 = u_A101;
        vk_A101 = v_A101;
        wk_A101 = w_A101;
        
        u_A000 = u_tilde_A000 + theta * div_p1_A000;
        v_A000 = v_tilde_A000 + theta * div_p2_A000;
        w_A000 = w_tilde_A000 + theta * div_p3_A000;
        
        u_A110 = u_tilde_A110 + theta * div_p1_A110;
        v_A110 = v_tilde_A110 + theta * div_p2_A110;
        w_A110 = w_tilde_A110 + theta * div_p3_A110;
        
        u_A011 = u_tilde_A011 + theta * div_p1_A011;
        v_A011 = v_tilde_A011 + theta * div_p2_A011;
        w_A011 = w_tilde_A011 + theta * div_p3_A011;
        
        u_A101 = u_tilde_A101 + theta * div_p1_A101;
        v_A101 = v_tilde_A101 + theta * div_p2_A101;
        w_A101 = w_tilde_A101 + theta * div_p3_A101;

        error = error + sum( (u_A000-uk_A000).*(u_A000-uk_A000) + ...
                             (v_A000-vk_A000).*(v_A000-vk_A000) + ...
                             (w_A000-wk_A000).*(w_A000-wk_A000) , 'all')/N;
        error = error + sum( (u_A110-uk_A110).*(u_A110-uk_A110) + ...
                             (v_A110-vk_A110).*(v_A110-vk_A110) + ...
                             (w_A110-wk_A110).*(w_A110-wk_A110) , 'all')/N;
        error = error + sum( (u_A011-uk_A011).*(u_A011-uk_A011) + ...
                             (v_A011-vk_A011).*(v_A011-vk_A011) + ...
                             (w_A011-wk_A011).*(w_A011-wk_A011), 'all')/N;
        error = error + sum( (u_A101-uk_A101).*(u_A101-uk_A101) + ...
                             (v_A101-vk_A101).*(v_A101-vk_A101) + ...
                             (w_A101-wk_A101).*(w_A101-wk_A101), 'all')/N;
                         
        [ux_A000, uy_A000, uz_A000] = gradient000hv(u_A000, u_A110, u_A011, u_A101, 0);
        [vx_A000, vy_A000, vz_A000] = gradient000hv(v_A000, v_A110, v_A011, v_A101, 0);
        [wx_A000, wy_A000, wz_A000] = gradient000hv(w_A000, w_A110, w_A011, w_A101, 0);
        
        [ux_A110, uy_A110, uz_A110] = gradient110hv(u_A110, u_A000, u_A011, u_A101, 0);
        [vx_A110, vy_A110, vz_A110] = gradient110hv(v_A110, v_A000, v_A011, v_A101, 0);
        [wx_A110, wy_A110, wz_A110] = gradient110hv(w_A110, w_A000, w_A011, w_A101, 0);
        
        [ux_A011, uy_A011, uz_A011] = gradient011hv(u_A011, u_A000, u_A110, u_A101, 0);
        [vx_A011, vy_A011, vz_A011] = gradient011hv(v_A011, v_A000, v_A110, v_A101, 0);
        [wx_A011, wy_A011, wz_A011] = gradient011hv(w_A011, w_A000, w_A110, w_A101, 0);
        
        [ux_A101, uy_A101, uz_A101] = gradient101hv(u_A101, u_A000, u_A110, u_A011, 0);
        [vx_A101, vy_A101, vz_A101] = gradient101hv(v_A101, v_A000, v_A110, v_A011, 0);
        [wx_A101, wy_A101, wz_A101] = gradient101hv(w_A101, w_A000, w_A110, w_A011, 0);
    
        taut = tau / theta;
        g1_A000 = sqrt(ux_A000 .* ux_A000 + uy_A000 .* uy_A000 + uz_A000 .* uz_A000);        
        g2_A000 = sqrt(vx_A000 .* vx_A000 + vy_A000 .* vy_A000 + vz_A000 .* vz_A000);             
        g3_A000 = sqrt(wx_A000 .* wx_A000 + wy_A000 .* wy_A000 + wz_A000 .* wz_A000);
        
        g1_A110 = sqrt(ux_A110 .* ux_A110 + uy_A110 .* uy_A110 + uz_A110 .* uz_A110);        
        g2_A110 = sqrt(vx_A110 .* vx_A110 + vy_A110 .* vy_A110 + vz_A110 .* vz_A110);             
        g3_A110 = sqrt(wx_A110 .* wx_A110 + wy_A110 .* wy_A110 + wz_A110 .* wz_A110);
        
        g1_A011 = sqrt(ux_A011 .* ux_A011 + uy_A011 .* uy_A011 + uz_A011 .* uz_A011);        
        g2_A011 = sqrt(vx_A011 .* vx_A011 + vy_A011 .* vy_A011 + vz_A011 .* vz_A011);             
        g3_A011 = sqrt(wx_A011 .* wx_A011 + wy_A011 .* wy_A011 + wz_A011 .* wz_A011);
        
        g1_A101 = sqrt(ux_A101 .* ux_A101 + uy_A101 .* uy_A101 + uz_A101 .* uz_A101);        
        g2_A101 = sqrt(vx_A101 .* vx_A101 + vy_A101 .* vy_A101 + vz_A101 .* vz_A101);             
        g3_A101 = sqrt(wx_A101 .* wx_A101 + wy_A101 .* wy_A101 + wz_A101 .* wz_A101);
        
        ng1_A000 = 1.0 + taut * (g1_A000 + g2_A000 + g3_A000);
        ng2_A000 = 1.0 + taut * (g1_A000 + g2_A000 + g3_A000);
        ng3_A000 = 1.0 + taut * (g1_A000 + g2_A000 + g3_A000);
        
        ng1_A110 = 1.0 + taut * (g1_A110 + g2_A110 + g3_A110);
        ng2_A110 = 1.0 + taut * (g1_A110 + g2_A110 + g3_A110);
        ng3_A110 = 1.0 + taut * (g1_A110 + g2_A110 + g3_A110);
        
        ng1_A011 = 1.0 + taut * (g1_A011 + g2_A011 + g3_A011);
        ng2_A011 = 1.0 + taut * (g1_A011 + g2_A011 + g3_A011);
        ng3_A011 = 1.0 + taut * (g1_A011 + g2_A011 + g3_A011);
        
        ng1_A101 = 1.0 + taut * (g1_A101 + g2_A101 + g3_A101);
        ng2_A101 = 1.0 + taut * (g1_A101 + g2_A101 + g3_A101);
        ng3_A101 = 1.0 + taut * (g1_A101 + g2_A101 + g3_A101);
        
        p1x_A000 = (p1x_A000 + taut * ux_A000) ./ ng1_A000;
        p1y_A000 = (p1y_A000 + taut * uy_A000) ./ ng1_A000;
        p1z_A000 = (p1z_A000 + taut * uz_A000) ./ ng1_A000;
        p2x_A000 = (p2x_A000 + taut * vx_A000) ./ ng2_A000;
        p2y_A000 = (p2y_A000 + taut * vy_A000) ./ ng2_A000;
        p2z_A000 = (p2z_A000 + taut * vz_A000) ./ ng2_A000;
        p3x_A000 = (p3x_A000 + taut * wx_A000) ./ ng3_A000;
        p3y_A000 = (p3y_A000 + taut * wy_A000) ./ ng3_A000;
        p3z_A000 = (p3z_A000 + taut * wz_A000) ./ ng3_A000;
        
        p1x_A110 = (p1x_A110 + taut * ux_A110) ./ ng1_A110;
        p1y_A110 = (p1y_A110 + taut * uy_A110) ./ ng1_A110;
        p1z_A110 = (p1z_A110 + taut * uz_A110) ./ ng1_A110;
        p2x_A110 = (p2x_A110 + taut * vx_A110) ./ ng2_A110;
        p2y_A110 = (p2y_A110 + taut * vy_A110) ./ ng2_A110;
        p2z_A110 = (p2z_A110 + taut * vz_A110) ./ ng2_A110;
        p3x_A110 = (p3x_A110 + taut * wx_A110) ./ ng3_A110;
        p3y_A110 = (p3y_A110 + taut * wy_A110) ./ ng3_A110;
        p3z_A110 = (p3z_A110 + taut * wz_A110) ./ ng3_A110;
        
        p1x_A011 = (p1x_A011 + taut * ux_A011) ./ ng1_A011;
        p1y_A011 = (p1y_A011 + taut * uy_A011) ./ ng1_A011;
        p1z_A011 = (p1z_A011 + taut * uz_A011) ./ ng1_A011;
        p2x_A011 = (p2x_A011 + taut * vx_A011) ./ ng2_A011;
        p2y_A011 = (p2y_A011 + taut * vy_A011) ./ ng2_A011;
        p2z_A011 = (p2z_A011 + taut * vz_A011) ./ ng2_A011;
        p3x_A011 = (p3x_A011 + taut * wx_A011) ./ ng3_A011;
        p3y_A011 = (p3y_A011 + taut * wy_A011) ./ ng3_A011;
        p3z_A011 = (p3z_A011 + taut * wz_A011) ./ ng3_A011;
        
        p1x_A101 = (p1x_A101 + taut * ux_A101) ./ ng1_A101;
        p1y_A101 = (p1y_A101 + taut * uy_A101) ./ ng1_A101;
        p1z_A101 = (p1z_A101 + taut * uz_A101) ./ ng1_A101;
        p2x_A101 = (p2x_A101 + taut * vx_A101) ./ ng2_A101;
        p2y_A101 = (p2y_A101 + taut * vy_A101) ./ ng2_A101;
        p2z_A101 = (p2z_A101 + taut * vz_A101) ./ ng2_A101;
        p3x_A101 = (p3x_A101 + taut * wx_A101) ./ ng3_A101;
        p3y_A101 = (p3y_A101 + taut * wy_A101) ./ ng3_A101;
        p3z_A101 = (p3z_A101 + taut * wz_A101) ./ ng3_A101;
    end
end

%clift = cmax;
clift = 0;

para.sigma = 1;

u_A000 = medfilt3(u_A000,[mf,mf,mf]);
v_A000 = medfilt3(v_A000,[mf,mf,mf]);
w_A000 = medfilt3(w_A000,[mf,mf,mf]);

u_A110 = medfilt3(u_A110,[mf,mf,mf]);
v_A110 = medfilt3(v_A110,[mf,mf,mf]);
w_A110 = medfilt3(w_A110,[mf,mf,mf]);

u_A011 = medfilt3(u_A011,[mf,mf,mf]);
v_A011 = medfilt3(v_A011,[mf,mf,mf]);
w_A011 = medfilt3(w_A011,[mf,mf,mf]);

u_A101 = medfilt3(u_A101,[mf,mf,mf]);
v_A101 = medfilt3(v_A101,[mf,mf,mf]);
w_A101 = medfilt3(w_A101,[mf,mf,mf]);

u = lift_hv_nodetail(u_A000, u_A110, u_A011, u_A101, clift, 'min');
v = lift_hv_nodetail(v_A000, v_A110, v_A011, v_A101, clift, 'min');
w = lift_hv_nodetail(w_A000, w_A110, w_A011, w_A101, clift, 'min');

p1x = lift_hv_nodetail(p1x_A000, p1x_A110, p1x_A011, p1x_A101, clift, 'min');
p1y = lift_hv_nodetail(p1y_A000, p1y_A110, p1y_A011, p1y_A101, clift, 'min');
p1z = lift_hv_nodetail(p1z_A000, p1z_A110, p1z_A011, p1z_A101, clift, 'min');

p2x = lift_hv_nodetail(p2x_A000, p2x_A110, p2x_A011, p2x_A101, clift, 'min');
p2y = lift_hv_nodetail(p2y_A000, p2y_A110, p2y_A011, p2y_A101, clift, 'min');
p2z = lift_hv_nodetail(p2z_A000, p2z_A110, p2z_A011, p2z_A101, clift, 'min');

p3x = lift_hv_nodetail(p3x_A000, p3x_A110, p3x_A011, p3x_A101, clift, 'min');
p3y = lift_hv_nodetail(p3y_A000, p3y_A110, p3y_A011, p3y_A101, clift, 'min');
p3z = lift_hv_nodetail(p3z_A000, p3z_A110, p3z_A011, p3z_A101, clift, 'min');

        
        
       