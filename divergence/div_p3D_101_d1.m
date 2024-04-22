function dp = div_p3D_101_d1(px_A000, px_A101, py_A000, py_A101, pz_A000, pz_A101)
%------------------------------------------------------------------------------
% Calculates divergence for color 101 in rectangular setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
[pxx_A101, ~, ~] = gradient101d1(px_A101, px_A000, cmin);
[~, pyy_A101, ~] = gradient101d1(py_A101, py_A000, cmin);
[~, ~, pzz_A101] = gradient101d1(pz_A101, pz_A000, cmin);

dp = pxx_A101 + pyy_A101 + pzz_A101;
