function [u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z] = ....
    dual_tvl1_of(I0, I1, ...
    u, v, w, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, ...
    tau, lambda, theta, warps, epsilon, max_it, verbose)
%------------------------------------------------------------------------------
% Calculates primal dual optical flow in 3D on classical voxel grids
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
I0 = double(I0);
I1 = double(I1);
[nx, ny, nz] = size(I0);
N = nx*ny*nz;
l_t = lambda * theta;

for i = 1:warps
    
    I1w = warpImage3D(I1, u, v, w);
    % compute the derivatives, smoothing considering the directions
    [Ix0, Iy0, Iz0] = centralDiffDir(I0);
    [Ix1, Iy1, Iz1] = centralDiffDir(I1);
    
    %warp the derivative according to the current displacement
    I1wx = warpImage3D(Ix1, u, v, w); 
    I1wy = warpImage3D(Iy1, u, v, w);
    I1wz = warpImage3D(Iz1, u, v, w);
    
    %estimate derivative direction by mixing i0 and i1 derivatives
    I1wx = 0.5 * (I1wx + Ix0); 
    I1wy = 0.5 * (I1wy + Iy0);
    I1wz = 0.5 * (I1wz + Iz0);
    
    Ix2 = I1wx.*I1wx;
    Iy2 = I1wy.*I1wy;
    Iz2 = I1wz.*I1wz;
    
    grad = Ix2 + Iy2 + Iz2;
    
    rho_c = I1w - I0 - I1wx.*u - I1wy.*v - I1wz.*w;
    
    n = 0; error = inf;
    
    u = medfilt3(u,[5,5,5]);
    v = medfilt3(v,[5,5,5]);
    w = medfilt3(w,[5,5,5]);
    
    while error > 0.1*epsilon && n < max_it
        n = n + 1;
        rho = rho_c + I1wx .* u + I1wy .* v  + I1wz .* w;
        rho_msk_1 = rho < -l_t * grad;
        d1 = l_t * rho_msk_1 .* I1wx;
        d2 = l_t * rho_msk_1 .* I1wy;
        d3 = l_t * rho_msk_1 .* I1wz;
        
        rho_msk_2 = rho > l_t * grad;
        d1 = d1 - l_t * rho_msk_2 .* I1wx;
        d2 = d2 - l_t * rho_msk_2 .* I1wy;
        d3 = d3 - l_t * rho_msk_2 .* I1wz;
        
        rho_msk_else = ones(size(rho_msk_1)) - rho_msk_1 - rho_msk_2;
        if min(rho_msk_else, [], 'all') >= 0
            d1 = d1 - rho_msk_else .* rho ./ grad .* I1wx;
            d2 = d2 - rho_msk_else .* rho ./ grad .* I1wy;
            d3 = d3 - rho_msk_else .* rho ./ grad .* I1wz;
            
            d1(isnan(d1)) = 0;
            d2(isnan(d2)) = 0;
            d3(isnan(d3)) = 0;
            
        else
            sprintf('something went wrong');
        end
        
        u_tilde = u + d1;
        v_tilde = v + d2;
        w_tilde = w + d3;
        
        div_p1 = div_p3D(p1x, p1y, p1z);
        div_p2 = div_p3D(p2x, p2y, p2z);
        div_p3 = div_p3D(p3x, p3y, p3z);
        
        error = 0.0;
        
        uk = u;
        vk = v;
        wk = w;

        u = u_tilde + theta * div_p1;
        v = v_tilde + theta * div_p2;
        w = w_tilde + theta * div_p3;
                
        error = error + sum( (u-uk).*(u-uk) + (v-vk).*(v-vk) + (w-wk).*(w-wk), 'all')/N;
        [ux, uy, uz] = imgradientxyz(u, 'intermediate');
        [vx, vy, vz] = imgradientxyz(v, 'intermediate');
        [wx, wy, wz] = imgradientxyz(w, 'intermediate');
        
        taut = tau / theta;
        g1 = sqrt(ux .* ux + uy .* uy + uz .* uz);        
        g2 = sqrt(vx .* vx + vy .* vy + vz .* vz);             
        g3 = sqrt(wx .* wx + wy .* wy + wz .* wz);
        
        ng1 = 1.0 + taut * (g1 + g2 + g3);
        ng2 = 1.0 + taut * (g1 + g2 + g3);
        ng3 = 1.0 + taut * (g1 + g2 + g3);
        
        p1x = (p1x + taut * ux) ./ ng1;
        p1y = (p1y + taut * uy) ./ ng1;
        p1z = (p1z + taut * uz) ./ ng1;
        p2x = (p2x + taut * vx) ./ ng2;
        p2y = (p2y + taut * vy) ./ ng2;
        p2z = (p2z + taut * vz) ./ ng2;
        p3x = (p3x + taut * wx) ./ ng3;
        p3y = (p3y + taut * wy) ./ ng3;
        p3z = (p3z + taut * wz) ./ ng3;
    end
end

u = medfilt3(u,[5,5,5]);
v = medfilt3(v,[5,5,5]);
w = medfilt3(w,[5,5,5]);
        
        
        