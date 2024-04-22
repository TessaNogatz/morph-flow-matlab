function [dp_A000, dp_A110, dp_A011, dp_A101] = div_p3D_hv(...
    px_A000, px_A110, px_A011, px_A101, ...
    py_A000, py_A110, py_A011, py_A101, ...
    pz_A000, pz_A110, pz_A011, pz_A101)
%------------------------------------------------------------------------------
% Calculates divergence for color 000, 110, 011 and 101 in dodecahedral setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
cmin = 0;

[pxx_A000, ~, ~] = gradient000hv(px_A000, px_A110, px_A011, px_A101,cmin);
[~, pyy_A000, ~] = gradient000hv(py_A000, py_A110, py_A011, py_A101, cmin);
[~, ~, pzz_A000] = gradient000hv(pz_A000, pz_A110, pz_A011, pz_A101, cmin);

dp_A000 = pxx_A000 + pyy_A000 + pzz_A000;

[pxx_A110, ~, ~] = gradient110hv(px_A110, px_A000, px_A011, px_A101, cmin);
[~, pyy_A110, ~] = gradient110hv(py_A110, py_A000, py_A011, py_A101, cmin);
[~, ~, pzz_A110] = gradient110hv(pz_A110, pz_A000, pz_A011, pz_A101, cmin);

dp_A110 = pxx_A110 + pyy_A110 + pzz_A110;

[pxx_A011, ~, ~] = gradient011hv(px_A011, px_A000, px_A110, px_A101,cmin);
[~, pyy_A011, ~] = gradient011hv(py_A011, py_A000, py_A110, py_A101, cmin);
[~, ~, pzz_A011] = gradient011hv(pz_A011, pz_A000, pz_A110, pz_A101, cmin);

dp_A011 = pxx_A011 + pyy_A011 + pzz_A011;

[pxx_A101, ~, ~] = gradient101hv(px_A101, px_A000, px_A110, px_A011, cmin);
[~, pyy_A101, ~] = gradient101hv(py_A101, py_A000, py_A110, py_A011, cmin);
[~, ~, pzz_A101] = gradient101hv(pz_A101, pz_A000, pz_A110, pz_A011, cmin);

dp_A101 = pxx_A101 + pyy_A101 + pzz_A101;