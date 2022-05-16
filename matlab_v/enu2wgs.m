function XYZ_WGS = enu2wgs(XYZ_ENU,LATP,LONP,HP)

[XWGSP,YWGSP,ZWGSP] = llh2xyz(LATP,LONP,HP);

% Transformation matrix from ENU to WGS
R = [matR(LONP+90,'z')*matR(90-LATP,'x'), [XWGSP; YWGSP; ZWGSP]; 0 0 0 1 ];

% Transformation from WGS to ENU
XYZ_WGS = R*XYZ_ENU;

end

function R = matR(a,par)
% Rotation matrix defined by Riegl, a is angle in degrees, par is char 
% defining around which axes to rotate

switch(par)
    case 'x'
        R = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
    case 'y'
        R = [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
    case 'z'
        R = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
    otherwise
        error('You must put parameter')
end

end