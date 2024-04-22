function [dx, dy, dz] = gradient011hv(A011, A000, A110, A101, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 011 in dodecahedral setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
[m011, n011, l011] = size(A011);
[m000, n000, l000] = size(A000);
[m110, n110, l110] = size(A110);
[m101, n101, l101] = size(A101);
 
cellVol = 16 / 9 * sqrt(2) ^ 2 * sqrt(3);
adjustedFaceArea = 2 / 3 ; 
% To the actual face area formula * 2 * sqrt(2) is missing. This is
% cancelled out with the factors from unit normal computation

% triangle and filled

if n011 == n000 % n even
    if l011 == l000 % m even
        sFeB = stripF3D(extB3D(A000, cmin));
        sLeR = stripL3D(extR3D(A000, cmin));
        sLeRsFeB = stripL3D(extR3D(sFeB, cmin));
        S1 = -A000 + sLeRsFeB;
        S2x = sLeR - sFeB;
        S2z = - S2x;
    elseif l011 == l000-1 % m odd
        sB = stripB3D(A000);
        sF = stripF3D(A000);
        sLeR = stripL3D(extR3D(A000, cmin));
        S1 = sB + stripF3D(sLeR);
        S2x = stripB3D(sLeR) - sF;
        S2z = -S2x;
    else
        error(' gradient011hv - sizes do not match ');
    end
elseif n011 == n000-1 % n odd
    if l011 == l000 % m even
        sR = stripR3D(A000);
        sL = stripL3D(A000);
        sFeB = stripF3D(extB3D(A000, cmin));
        S1 = -sR + stripL3D(sFeB);
        S2x = +sL - stripR3D(sFeB);
        S2z = -S2x;
    elseif l011 == l000-1 % m odd
        sF = stripF3D(A000);
        sB = stripB3D(A000);
        sRsB = stripR3D(sB);
        sLsF = stripL3D(sF);
        sRsF = stripR3D(sF);
        sLsB = stripL3D(sB);
        S1 = -sRsB + sLsF;
        S2x = sLsB - sRsF;
        S2z = -S2x;
    else
        error(' gradient011hv - sizes do not match ');
    end
else 
    error(' gradient011hv - sizes do not match ');
end

% triangle with rhomboid

if m011 == m110 % m even
    if l011 == l110 % l even
        sFeB = stripF3D(extB3D(A110, cmin));
        sDeU = stripD3D(extU3D(A110, cmin));
        sDeUsFeB = stripD3D(extU3D(sFeB, cmin));
        T1y = A110 - sDeUsFeB;
        T1z = -T1y;
        T2 = sFeB - sDeU;
    elseif l011 == l110-1 % l odd
       sB = stripB3D(A110);
        sF = stripF3D(A110);
        sDeU = stripD3D(extU3D(A110, cmin));
        T1y = sB - stripF3D(sDeU);
        T1z = -T1y;
        T2 = sF - stripB3D(sDeU);
    else
        error(' gradient011hv - sizes do not match ');
    end
elseif m011 == m110+1 % m odd
    if l011 == l110 % l even
        eD = extD3D(A110, cmin);
        eU = extU3D(A110, cmin);
        sFeB = stripF3D(extB3D(A110, cmin));
        T1y = eD - extU3D(sFeB, cmin);
        T1z = -T1y;
        T2 = extD3D(sFeB, cmin) - eU;
    elseif l011 == l110-1 % l odd
        eU = extU3D(A110, cmin);
        eD = extD3D(A110, cmin);
        eDsB = stripB3D(eD);
        eUsF = stripF3D(eU);
        sFeD = stripF3D(eD);
        sBeU = stripB3D(eU);
        T1y = eDsB - eUsF;
        T1z = -T1y;
        T2 = sFeD - sBeU;
    else
        error(' gradient011hv - sizes do not match ');
    end
else 
    error(' gradient011hv - sizes do not match ');
end

% triangle with star

if n011 == n101 % n even
    if m011 == m101 % m even
        sLeR = stripL3D(extR3D(A101, cmin));
        sDeU = stripD3D(extU3D(A101, cmin));
        sLeRsDeU = stripL3D(extR3D(sDeU, cmin));
        R1x = -A101 + sLeRsDeU;
        R1y = -R1x;
        R2  = sLeR-sDeU;
    elseif m011 == m101+1 % m odd
        eD = extD3D(A101, cmin);
        eU = extU3D(A101, cmin);
        sLeR = stripL3D(extR3D(A101, cmin));
        R1x = -eD + extU3D(sLeR, cmin);
        R1y = -R1x;
        R2 = extD3D(sLeR, cmin) - eU; 
    else
        error(' gradient011hv - sizes do not match ');
    end
elseif n011 == n101-1 % n odd
    if m011 == m101 % m even
        sR = stripR3D(A101);
        sL = stripL3D(A101);
        sDeU = stripD3D(extU3D(A101, cmin));
        R1x = -sR + stripL3D(sDeU);
        R1y = -R1x;
        R2  = sL -stripR3D(sDeU);
    elseif m011 == m101+1 % m odd
        eD = extD3D(A101, cmin);
        eU = extU3D(A101, cmin);
        sReD = stripR3D(eD);
        sLeU = stripL3D(eU);
        sLeD = stripL3D(eD);
        sReU = stripR3D(eU);
        R1x = -sReD + sLeU;
        R1y = -R1x;
        R2 = sLeD - sReU;
    else
        error(' gradient011hv - sizes do not match ');
    end
else 
    error(' gradient011hv - sizes do not match ');
end

dx = adjustedFaceArea / cellVol * ( S1 + S2x + R1x + R2);
dy = adjustedFaceArea / cellVol * ( R1y + R2 + T1y + T2);
dz = adjustedFaceArea / cellVol * ( S1 + S2z + T1z + T2);
