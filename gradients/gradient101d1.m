function [grad101x, grad101y, grad101z] = gradient101d1(A101, A000, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 101 in rectangular setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
[m101, ~, l101] = size(A101);
[m000, ~, l000] = size(A000);

cellVol = sqrt(2) * sqrt(2); 

grad101x = (1/cellVol) * 2 * 0.5 * ... % face area is sqrt 2 * sqrt 2
    (stripL3D(extR3D(A101, cmin)) - stripR3D(extL3D(A101, cmin)));

if l101 == l000 % n even
    if m101 == m000 % l even
        FUBD = - A000 + extB3D(stripF3D(extD3D(stripU3D(A000), cmin)), cmin);
        FDBUz = - stripU3D(extD3D(A000, cmin)) + stripF3D(extB3D(A000, cmin));
        FDBUy = + stripU3D(extD3D(A000, cmin)) - stripF3D(extB3D(A000, cmin));
    elseif m101 == m000-1 % l odd
        FUBD = - stripD3D(A000) + extB3D(stripF3D(stripU3D(A000)), cmin);
        FDBUz = - stripU3D(A000) + stripD3D(stripF3D(extB3D(A000, cmin)));
        FDBUy = + stripU3D(A000) - stripD3D(stripF3D(extB3D(A000, cmin)));
    else
        error(' d2Lift101max - sizes do not match ');
    end
elseif l101 == l000-1 % l odd
    if m101 == m000 % m even
        FUBD = - stripB3D(A000) + stripF3D(extD3D(stripU3D(A000), cmin));
        FDBUz = - stripB3D(stripU3D(extD3D(A000, cmin))) + stripF3D(A000);
        FDBUy = + stripB3D(stripU3D(extD3D(A000, cmin))) - stripF3D(A000);
    elseif m101 == m000-1 % l odd
        FUBD = - stripB3D(stripD3D(A000)) + stripF3D(stripU3D(A000));
        FDBUy = + stripB3D(stripU3D(A000)) - stripD3D(stripF3D(A000));
        FDBUz = - stripB3D(stripU3D(A000)) + stripD3D(stripF3D(A000));
    else
        error(' d2Lift101max - sizes do not match ');
    end
else 
    error(' d2Lift101max - sizes do not match ');
end

grad101y = (1/cellVol) * 2 * sqrt(2) * 0.5 * ( FUBD + FDBUy );
grad101z = (1/cellVol) * 2 * sqrt(2) * 0.5 * ( FUBD + FDBUz );
