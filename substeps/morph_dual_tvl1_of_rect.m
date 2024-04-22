function [u_A000, u_A101, v_A000, v_A101, w_A000, w_A101, ...
    p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ...
    morph_dual_tvl1_of_rect(I0_A000, I0_A101, I1_A000, I1_A101, ...
    u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, ...
    tau, lambda, theta, warps, epsilon, max_it, verbose)

%------------------------------------------------------------------------------
% Calculates primal dual optical flow in 3D on rectangular voxel grids
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
[nx, ny, nz] = size(I0_A000);
N = nx*2*ny*2*nz*2;
l_t = lambda * theta;

minO = min(min(u, [], 'all'), 0);
maxO = max(max(u, [], 'all'), 0);
cmin = minO-(maxO-minO);
cmax = maxO+(maxO-minO);

yA101 = size(I0_A101,1);
zA101 = size(I0_A101,3);

[u_A000, u_A101] = lift_d2_nodetail(u, yA101, zA101, cmax, 'min');
[v_A000, v_A101] = lift_d2_nodetail(v, yA101, zA101, cmax, 'min');
[w_A000, w_A101] = lift_d2_nodetail(w, yA101, zA101, cmax, 'min');

[p1x_A000, p1x_A101] = lift_d2_nodetail(p1x, yA101, zA101, cmax, 'min');
[p1y_A000, p1y_A101] = lift_d2_nodetail(p1y, yA101, zA101, cmax, 'min');
[p1z_A000, p1z_A101] = lift_d2_nodetail(p1z, yA101, zA101, cmax, 'min');

[p2x_A000, p2x_A101] = lift_d2_nodetail(p2x, yA101, zA101, cmax, 'min');
[p2y_A000, p2y_A101] = lift_d2_nodetail(p2y, yA101, zA101, cmax, 'min');
[p2z_A000, p2z_A101] = lift_d2_nodetail(p2z, yA101, zA101, cmax, 'min');

[p3x_A000, p3x_A101] = lift_d2_nodetail(p3x, yA101, zA101, cmax, 'min');
[p3y_A000, p3y_A101] = lift_d2_nodetail(p3y, yA101, zA101, cmax, 'min');
[p3z_A000, p3z_A101] = lift_d2_nodetail(p3z, yA101, zA101, cmax, 'min');

for i = 1:warps
    
    I1w_A000 = warpImage3D(I1_A000, u_A000, v_A000, w_A000);
    I1w_A101 = warpImage3D(I1_A101, u_A101, v_A101, w_A101);
    % compute the derivatives, smoothing considering the directions
    [Ix0_A000, Iy0_A000, Iz0_A000] = gradient000d1(I0_A000, I0_A101, 0);
    [Ix0_A101, Iy0_A101, Iz0_A101] = gradient101d1(I0_A101, I0_A000, 0);
    [Ix1_A000, Iy1_A000, Iz1_A000] = gradient000d1(I1_A000, I1_A101, 0);
    [Ix1_A101, Iy1_A101, Iz1_A101] = gradient101d1(I1_A101, I1_A000, 0);
    %[Ix1, Iy1, Iz1] = centralDiffDir(I1);
    
    %warp the derivative according to the current displacement
    I1wx_A000 = warpImage3D(Ix1_A000, u_A000, v_A000, w_A000); 
    I1wx_A101 = warpImage3D(Ix1_A101, u_A101, v_A101, w_A101);
    I1wy_A000 = warpImage3D(Iy1_A000, u_A000, v_A000, w_A000); 
    I1wy_A101 = warpImage3D(Iy1_A101, u_A101, v_A101, w_A101);
    I1wz_A000 = warpImage3D(Iz1_A000, u_A000, v_A000, w_A000); 
    I1wz_A101 = warpImage3D(Iz1_A101, u_A101, v_A101, w_A101);
     
    %estimate derivative direction by mixing i0 and i1 derivatives
    I1wx_A000 = 0.5 * (I1wx_A000 + Ix0_A000); 
    I1wx_A101 = 0.5 * (I1wx_A101 + Ix0_A101); 
    I1wy_A000 = 0.5 * (I1wy_A000 + Iy0_A000); 
    I1wy_A101 = 0.5 * (I1wy_A101 + Iy0_A101); 
    I1wz_A000 = 0.5 * (I1wz_A000 + Iz0_A000); 
    I1wz_A101 = 0.5 * (I1wz_A101 + Iz0_A101); 
    
    Ix2_A000 = I1wx_A000.*I1wx_A000;
    Ix2_A101 = I1wx_A101.*I1wx_A101;
    Iy2_A000 = I1wy_A000.*I1wy_A000;
    Iy2_A101 = I1wy_A101.*I1wy_A101;
    Iz2_A000 = I1wz_A000.*I1wz_A000;
    Iz2_A101 = I1wz_A101.*I1wz_A101;
    
    grad_A000 = Ix2_A000 + Iy2_A000 + Iz2_A000;
    grad_A101 = Ix2_A101 + Iy2_A101 + Iz2_A101;
    
    rho_c_A000 = I1w_A000 - I0_A000 - ...
        I1wx_A000.*u_A000 - I1wy_A000.*v_A000 - I1wz_A000.*w_A000;
    rho_c_A101 = I1w_A101 - I0_A101 - ...
        I1wx_A101.*u_A101 - I1wy_A101.*v_A101 - I1wz_A101.*w_A101;
    
    n = 0; error = inf;
       
    mf = 3;

    u_A000 = medfilt3(u_A000,[mf, mf, mf]);
    v_A000 = medfilt3(v_A000,[mf, mf, mf]);
    w_A000 = medfilt3(w_A000,[mf, mf, mf]);
        
    u_A101 = medfilt3(u_A101,[mf, mf, mf]);
    v_A101 = medfilt3(v_A101,[mf, mf, mf]);
    w_A101 = medfilt3(w_A101,[mf, mf, mf]);
       

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
        
        u_tilde_A101 = u_A101 + d1_A101;
        v_tilde_A101 = v_A101 + d2_A101;
        w_tilde_A101 = w_A101 + d3_A101;
        
        [div_p1_A000, div_p1_A101] = div_p3D_d1(...
            p1x_A000, p1x_A101, p1y_A000, p1y_A101, p1z_A000, p1z_A101);
        [div_p2_A000, div_p2_A101] = div_p3D_d1(...
            p2x_A000, p2x_A101, p2y_A000, p2y_A101, p2z_A000, p2z_A101);
        [div_p3_A000, div_p3_A101] = div_p3D_d1(...
            p3x_A000, p3x_A101, p3y_A000, p3y_A101, p3z_A000, p3z_A101);
        
        error = 0.0;
        
        uk_A000 = u_A000;
        vk_A000 = v_A000;
        wk_A000 = w_A000;
        
        uk_A101 = u_A101;
        vk_A101 = v_A101;
        wk_A101 = w_A101;
        
        u_A000 = u_tilde_A000 + theta * div_p1_A000;
        v_A000 = v_tilde_A000 + theta * div_p2_A000;
        w_A000 = w_tilde_A000 + theta * div_p3_A000;
        
        u_A101 = u_tilde_A101 + theta * div_p1_A101;
        v_A101 = v_tilde_A101 + theta * div_p2_A101;
        w_A101 = w_tilde_A101 + theta * div_p3_A101;
        
        error = error + sum( (u_A000-uk_A000).*(u_A000-uk_A000) + ...
                             (v_A000-vk_A000).*(v_A000-vk_A000) + ...
                             (w_A000-wk_A000).*(w_A000-wk_A000), 'all')/N; 
        error = error + sum( (u_A101-uk_A101).*(u_A101-uk_A101) + ...
                             (v_A101-vk_A101).*(v_A101-vk_A101) + ...
                             (w_A101-wk_A101).*(w_A101-wk_A101), 'all')/N;
                         
        
        
        [ux_A000, uy_A000, uz_A000] = gradient000d1(u_A000, u_A101, 0);
        [vx_A000, vy_A000, vz_A000] = gradient000d1(v_A000, v_A101, 0);
        [wx_A000, wy_A000, wz_A000] = gradient000d1(w_A000, w_A101, 0);
        
        [ux_A101, uy_A101, uz_A101] = gradient101d1(u_A101, u_A000, 0);
        [vx_A101, vy_A101, vz_A101] = gradient101d1(v_A101, v_A000, 0);
        [wx_A101, wy_A101, wz_A101] = gradient101d1(w_A101, w_A000, 0);
        
        taut = tau / theta;
        g1_A000 = sqrt(ux_A000 .* ux_A000 + uy_A000 .* uy_A000 + ...
            uz_A000 .* uz_A000);        
        g2_A000 = sqrt(vx_A000 .* vx_A000 + vy_A000 .* vy_A000 + ...
            vz_A000 .* vz_A000);             
        g3_A000 = sqrt(wx_A000 .* wx_A000 + wy_A000 .* wy_A000 + ...
            wz_A000 .* wz_A000);
        
        g1_A101 = sqrt(ux_A101 .* ux_A101 + uy_A101 .* uy_A101 + ...
            uz_A101 .* uz_A101);        
        g2_A101 = sqrt(vx_A101 .* vx_A101 + vy_A101 .* vy_A101 + ...
            vz_A101 .* vz_A101);             
        g3_A101 = sqrt(wx_A101 .* wx_A101 + wy_A101 .* wy_A101 + ...
            wz_A101 .* wz_A101);
        
        ng1_A000 = 1.0 + taut * (g1_A000 + g2_A000 + g3_A000);
        ng2_A000 = 1.0 + taut * (g2_A000 + g2_A000 + g3_A000);
        ng3_A000 = 1.0 + taut * (g3_A000 + g2_A000 + g3_A000);
        
        ng1_A101 = 1.0 + taut * (g1_A101 + g2_A101 + g3_A101);
        ng2_A101 = 1.0 + taut * (g2_A101 + g2_A101 + g3_A101);
        ng3_A101 = 1.0 + taut * (g3_A101 + g2_A101 + g3_A101);
        
        p1x_A000 = (p1x_A000 + taut * ux_A000) ./ ng1_A000;
        p1y_A000 = (p1y_A000 + taut * uy_A000) ./ ng1_A000;
        p1z_A000 = (p1z_A000 + taut * uz_A000) ./ ng1_A000;
        p2x_A000 = (p2x_A000 + taut * vx_A000) ./ ng2_A000;
        p2y_A000 = (p2y_A000 + taut * vy_A000) ./ ng2_A000;
        p2z_A000 = (p2z_A000 + taut * vz_A000) ./ ng2_A000;
        p3x_A000 = (p3x_A000 + taut * wx_A000) ./ ng3_A000;
        p3y_A000 = (p3y_A000 + taut * wy_A000) ./ ng3_A000;
        p3z_A000 = (p3z_A000 + taut * wz_A000) ./ ng3_A000;
        
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

p1x = {p1x_A000, p1x_A101};
p1y = {p1y_A000, p1y_A101};
p1z = {p1z_A000, p1z_A101};

p2x = {p2x_A000, p2x_A101};
p2y = {p2y_A000, p2y_A101};
p2z = {p2z_A000, p2z_A101};

p3x = {p3x_A000, p3x_A101};
p3y = {p3y_A000, p3y_A101};
p3z = {p3z_A000, p3z_A101};

        
 u_A000 = medfilt3(u_A000,[mf, mf, mf]);
 v_A000 = medfilt3(v_A000,[mf, mf, mf]);
 w_A000 = medfilt3(w_A000,[mf, mf, mf]);
        
 u_A101 = medfilt3(u_A101,[mf, mf, mf]);
 v_A101 = medfilt3(v_A101,[mf, mf, mf]);
 w_A101 = medfilt3(w_A101,[mf, mf, mf]);
       
        
       