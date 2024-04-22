function dp = div_p3D(px, py, pz)
%------------------------------------------------------------------------------
% Calculates divergence of a 3D vectorfield
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
pxx = zeros(size(px)); pyy = zeros(size(py)); pzz = zeros(size(pz));

pyy(2:end, :,:) = diff(py, 1, 1);
pxx(:,2:end,:) = diff(px, 1, 2);
pzz(:,:,2:end) = diff(pz, 1, 3);

dp = pxx + pyy + pzz;

dp(2:end,1,1) = px(2:end,1,1) + py(2:end,1,1) - py(1:end-1,1,1) ...
                              + pz(2:end,1,1) - pz(1:end-1,1,1);
dp(1,2:end,1) = py(1,2:end,1) + px(1,2:end,1) - px(1,1:end-1,1)...
                              + pz(1,2:end,1) - pz(1,1:end-1,1);
dp(1,1,2:end) = pz(1,1,2:end) + px(1,1,2:end) - px(1,1,1:end-1)...
                              + py(1,1,2:end) - py(1,1,1:end-1);

dp(1,1) = px(1,1) + py(1,1) + pz(1,1);