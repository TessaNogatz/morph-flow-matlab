function dp = div_p3D_000_d1(px_A000, px_A101, py_A000, py_A101, pz_A000, pz_A101)
%------------------------------------------------------------------------------
% Calculates divergence for color 000 in rectangular setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
[pxx_A000, ~, ~] = gradient000d1(px_A000, px_A101, cmin);
[~, pyy_A000, ~] = gradient000d1(py_A000, py_A101, cmin);
[~, ~, pzz_A000] = gradient000d1(pz_A000, pz_A101, cmin);

dp = pxx_A000 + pyy_A000 + pzz_A000;
