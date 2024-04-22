function [dx, dy, dz] = gradient101hv(A101, A000, A110, A011, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 101 in dodecahedral setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 

[m101, n101, l101] = size(A101);
[m000, n000, l000] = size(A000);
[m110, n110, l110] = size(A110);
[m011, n011, l011] = size(A011);

 
cellVol = 16 / 9 * sqrt(2) ^ 2 * sqrt(3);
adjustedFaceArea = 2 / 3 ; 
% To the actual face area formula * 2 * sqrt(2) is missing. This is
% cancelled out with the factors from unit normal computation

% star and filled

if l101 == l000 % n even
        if m101 == m000 % l even
            sUeD = stripU3D(extD3D(A000, cmin));
            sFeB = stripF3D(extB3D(A000, cmin));
            eBsFeDsU = extB3D(stripF3D(sUeD), cmin);
            S1 = -A000 + eBsFeDsU;
            S2z = -sUeD + sFeB;
            S2y = -S2z;
        elseif m101 == m000-1 % l odd
            sD = stripD3D(A000);
            sU = stripU3D(A000);
            sFeB = stripF3D(extB3D(A000, cmin));
            S1 = -sD + stripU3D(sFeB);
            S2z = -sU + stripD3D(sFeB);
            S2y = -S2z;
        else
            error(' d2Lift101max - sizes do not match ');
        end
    elseif l101 == l000-1 % l odd
        if m101 == m000 % m even
            sB = stripB3D(A000);
            sF = stripF3D(A000);
            eDsU = extD3D(stripU3D(A000), cmin);
            S1 = -sB + stripF3D(eDsU);
            S2z = -stripB3D(eDsU) + sF;
            S2y = -S2z;
        elseif m101 == m000-1 % l odd
            sD = stripD3D(A000);
            sU = stripU3D(A000);
            sBsD = stripB3D(sD);
            sFsU = stripF3D(sU);
            sBsU = stripB3D(sU);
            S1 = -sBsD + sFsU;
            S2z = -sBsU + sBsD;
            S2y = -S2z;
        else
            error(' d2Lift101max - sizes do not match ');
        end
    else 
        error(' d2Lift101max - sizes do not match ');
end

    
% star and rhomboid

if n101 == n110 % n even
    if l101 == l110 % l even
        sReL = stripR3D(extL3D(A110, cmin));
        sFeB = stripF3D(extB3D(A110, cmin));
        sFeBsReL = stripF3D(extB3D(sReL, cmin));
        R1x = -sFeBsReL + A110;
        R1z = -R1x;
        R2 = sReL + sFeB;
    elseif l101 == l110-1 % l odd
        sB = stripB3D(A110);
        sF = stripF3D(A110);
        sReL = stripR3D(extL3D(A110, cmin));
        R1x = -stripF3D(sReL) + sB;
        R1z = -R1x;
        R2 = - stripB3D(sReL) + sF;
    else
        error(' d1Lift101max - sizes do not match ');
    end
elseif n101 == n110+1 % n odd
    if l101 == l110 % l even
        eR = extR3D(A110, cmin);
        eL = extL3D(A110, cmin);
        sFeB = stripF3D(extB3D(A110, cmin));
        R1x = -extL3D(sFeB, cmin) + eR;
        R1z = -R1x;
        R2 = -eL + extR3D(sFeB, cmin);
    elseif l101 == l110-1 % l odd
       eL = extL3D(A110, cmin);
        eR = extR3D(A110, cmin);
        sFeL = stripF3D(eL);
        sBeR = stripB3D(eR);
        eLsB = stripB3D(eL);
        eRsF = stripF3D(eR);
        R1x = -sFeL + sBeR;
        R1z = -R1x;
        R2 = -eLsB + eRsF;
    else
        error(' d1Lift101max - sizes do not match ');
    end
else 
    error(' d1Lift101max - sizes do not match ');
end

% star and triangle

if n101 == n011 % n even
    if m101 == m011 % m even
        sUeD = stripU3D(extD3D(A011, cmin));
        sReL = stripR3D(extL3D(A011, cmin));
        eLsReDsU = stripR3D(extL3D(sUeD, cmin));
        T1x = A011 - eLsReDsU;
        T1y = -T1x;
        T2 = sUeD - sReL;
    elseif m101 == m011-1 % m odd
        sD =  stripD3D(A011);
        sU = stripU3D(A011);
        eLsR = stripR3D(extL3D(A011, cmin));
        eLsRsU = stripU3D(eLsR);
        sDsReL = stripD3D(eLsR);
        T1x = sD - eLsRsU;
        T1y = -T1x;
        T2 = sU - sDsReL;
    else
        error(' d1Lift101max - sizes do not match ');
    end
elseif n101 == n011+1 % n odd
    if m101 == m011 % m even
        eR = extR3D(A011, cmin);
        eL = extL3D(A011, cmin);
        eDsU = extD3D(stripU3D(A011), cmin);
        eLeDsU = extL3D(eDsU, cmin);
        eReDsU = extR3D(eDsU, cmin);
        T1x = eR - eLeDsU;
        T1y = -T1x;
        T2 = eReDsU - eL;
    elseif m101 == m011-1 % m odd
        sU = stripU3D(A011);
        sD = stripD3D(A011);
        eRsD = extR3D(sD, cmin);
        eLsU = extL3D(sU, cmin);
        eRsU = extR3D(sU, cmin);
        eLsD = extL3D(sD, cmin);
        T1x = eRsD - eLsU;
        T1y = -T1x;
        T2 = eRsU - eLsD;
    else
        error(' d1Lift101max - sizes do not match ');
    end
else 
    error(' d1Lift101max - sizes do not match ');
end

dx = adjustedFaceArea / cellVol * ( R1x + R2 + T1x + T2);
dy = adjustedFaceArea / cellVol * ( S1 + S2y + T1y + T2);
dz = adjustedFaceArea / cellVol * ( S1 + S2z + R1z + R2);
