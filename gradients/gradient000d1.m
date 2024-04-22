function [grad000x, grad000y, grad000z] = gradient000d1(A000, A101, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 000 in rectangular setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
[m000, ~, l000] = size(A000);
[m101, ~, l101] = size(A101);

cellVol = sqrt(2) * sqrt(2); 

grad000x = (1/cellVol) * 2 * 0.5 * ... % face area is sqrt 2 * sqrt 2
    (stripL3D(extR3D(A000, cmin)) - stripR3D(extL3D(A000, cmin)));

if l000 == l101 % n even
    if m000 == m101 % m even
        BUFDz = - extF3D(stripB3D(A101), cmin) + stripD3D(extU3D(A101, cmin));
        BUFDy =   extF3D(stripB3D(A101), cmin) - stripD3D(extU3D(A101, cmin));
        BDFU = A101 - extF3D(stripB3D(extU3D(stripD3D(A101), cmin)), cmin);
    elseif m000 == m101+1 % m odd
        BUFDz = - extF3D(stripB3D(extD3D(A101, cmin)), cmin) + extU3D(A101, cmin);
        BUFDy =   extF3D(stripB3D(extD3D(A101, cmin)), cmin) - extU3D(A101, cmin);
        BDFU = extD3D(A101, cmin) - extF3D(stripB3D(extU3D(A101, cmin)), cmin);
    else
        error(' d1Lift000max - sizes do not match ');
    end
elseif l000 == l101+1 % n odd
    if m000 == m101 % m even
        BUFDz = - extF3D(A101, cmin) + extB3D(extU3D(stripD3D(A101), cmin), cmin);
        BUFDy = + extF3D(A101, cmin) - extB3D(extU3D(stripD3D(A101), cmin), cmin);
        BDFU = extB3D(A101, cmin) - extF3D(stripD3D(extU3D(A101, cmin)), cmin);
    elseif m000 == m101+1 % m odd
        BUFDz = - extF3D(extD3D(A101, cmin), cmin) + extB3D(extU3D(A101, cmin), cmin);
        BUFDy =   extF3D(extD3D(A101, cmin), cmin) - extB3D(extU3D(A101, cmin), cmin);
        BDFU = extB3D(extD3D(A101, cmin), cmin) - extF3D(extU3D(A101, cmin), cmin);
    else
        error(' d1Lift000max - sizes do not match ');
    end
else 
    error(' d1Lift000max - sizes do not match ');
end

grad000y = (1/cellVol) * 2 * sqrt(2) * 0.5 * ( BDFU + BUFDy );
grad000z = (1/cellVol) * 2 * sqrt(2) * 0.5 * ( BDFU + BUFDz );

