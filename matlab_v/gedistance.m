function [s12, azi1, azi2, S12] = gedistance(lat1, lon1, lat2, lon2, ellipsoid)
%GEDISTANCE  Great ellipse distance on an ellipsoid
%
%   [s12, azi1, azi2] = GEDISTANCE(lat1, lon1, lat2, lon2)
%   [s12, azi1, azi2, S12] = GEDISTANCE(lat1, lon1, lat2, lon2, ellipsoid)
%
%   solves the inverse great ellipse problem of finding of length and
%   azimuths of the great ellipse between points specified by lat1, lon1,
%   lat2, lon2.  The input latitudes and longitudes, lat1, lon1, lat2,
%   lon2, can be scalars or arrays of equal size and must be expressed in
%   degrees.  The ellipsoid vector is of the form [a, e], where a is the
%   equatorial radius in meters, e is the eccentricity.  If ellipsoid is
%   omitted, the WGS84 ellipsoid (more precisely, the value returned by
%   defaultellipsoid) is used.  The output s12 is the distance in meters
%   and azi1 and azi2 are the forward azimuths at the end points in
%   degrees.  The optional output S12 is the area between the great ellipse
%   and the equator (in meters^2).  gedoc gives an example and provides
%   additional background information.  gedoc also gives the restrictions
%   on the allowed ranges of the arguments.
%
%   When given a combination of scalar and array inputs, the scalar inputs
%   are automatically expanded to match the size of the arrays.
%
%   geoddistance solves the equivalent geodesic problem and usually this is
%   preferable to using GEDISTANCE.
%
%   See also GEDOC, GERECKON, DEFAULTELLIPSOID, GEODDISTANCE, GEODRECKON.

% Copyright (c) Charles Karney (2014-2015) <charles@karney.com>.
%
% This file was distributed with GeographicLib 1.42.

  narginchk(4, 5)
  if nargin < 5, ellipsoid = defaultellipsoid; end
  try
    S = size(lat1 + lon1 + lat2 + lon2);
  catch
    error('lat1, lon1, s12, azi1 have incompatible sizes')
  end
  if length(ellipsoid(:)) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  Z = zeros(S);
  lat1 = lat1 + Z; lon1 = lon1 + Z;
  lat2 = lat2 + Z; lon2 = lon2 + Z;

  degree = pi/180;
  tiny = sqrt(realmin);

  a = ellipsoid(1);
  e2 = ellipsoid(2)^2;
  f = e2 / (1 + sqrt(1 - e2));

  f1 = 1 - f;

  areap = nargout >= 4;

  lon12 = AngDiff(AngNormalize(lon1(:)), AngNormalize(lon2(:)));
  lon12 = AngRound(lon12);

  phi = lat1 * degree;
  sbet1 = f1 * sin(phi); cbet1 = cos(phi); cbet1(lat1 == -90) = tiny;
  [sbet1, cbet1] = norm2(sbet1, cbet1);

  phi = lat2 * degree;
  sbet2 = f1 * sin(phi); cbet2 = cos(phi); cbet2(abs(lat2) == 90) = tiny;
  [sbet2, cbet2] = norm2(sbet2, cbet2);

  lam12 = lon12 * degree;
  slam12 = sin(lam12); slam12(lon12 == 180) = 0; clam12 = cos(lam12);

  % Solve great circle
  sgam1 = cbet2 .* slam12; cgam1 = +cbet1 .* sbet2 - sbet1 .* cbet2 .* clam12;
  sgam2 = cbet1 .* slam12; cgam2 = -sbet1 .* cbet2 + cbet1 .* sbet2 .* clam12;
  ssig12 = hypot(sgam1, cgam1);
  csig12 = sbet1 .* sbet2 + cbet1 .* cbet2 .* clam12;
  [sgam1, cgam1] = norm2(sgam1, cgam1);
  [sgam2, cgam2] = norm2(sgam2, cgam2);
  % no need to normalize [ssig12, csig12]

  cgam0 = hypot(cgam1, sgam1 .* sbet1);

  ssig1 = sbet1; csig1 = cbet1 .* cgam1;
  [ssig1, csig1] = norm2(ssig1, csig1);
  ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
  csig2 = csig1 .* csig12 - ssig1 .* ssig12;

  k2 = e2 * cgam0.^2;
  epsi = k2 ./ (2 * (1 + sqrt(1 - k2)) - k2);
  C1a = C1f(epsi);
  A1 = a * (1 + A1m1f(epsi)) .* (1 - epsi)./(1 + epsi);
  s12 = A1 .* (atan2(ssig12, csig12) + ...
               (SinCosSeries(true, ssig2, csig2, C1a) - ...
                SinCosSeries(true, ssig1, csig1, C1a)));
  azi1 = atan2(sgam1, cgam1 .* sqrt(1 - e2 * cbet1.^2)) / degree;
  azi2 = atan2(sgam2, cgam2 .* sqrt(1 - e2 * cbet2.^2)) / degree;

  s12 = reshape(s12, S); azi1 = reshape(azi1, S); azi2 = reshape(azi2, S);

  if areap
    sgam0 = sgam1 .* cbet1;
    A4 = (a^2 * e2) * cgam0 .* sgam0;

    n = f / (2 - f);
    G4x = G4coeff(n);
    G4a = C4f(epsi, G4x);
    B41 = SinCosSeries(false, ssig1, csig1, G4a);
    B42 = SinCosSeries(false, ssig2, csig2, G4a);
    S12 = A4 .* (B42 - B41);
    S12(cgam0 == 0 | sgam0 == 0) = 0;

    sgam12 = sgam2 .* cgam1 - cgam2 .* sgam1;
    cgam12 = cgam2 .* cgam1 + sgam2 .* sgam1;
    s = sgam12 == 0 & cgam12 < 0;
    sgam12(s) = tiny * cgam1(s); cgam12(s) = -1;
    gam12 = atan2(sgam12, cgam12);

    l = abs(gam12) < 1;
    dlam12 = 1 + clam12(l); dbet1 = 1 + cbet1(l); dbet2 = 1 + cbet2(l);
    gam12(l) = ...
        2 * atan2(slam12(l) .* (sbet1(l) .* dbet2 + sbet2(l) .* dbet1), ...
                  dlam12    .* (sbet1(l) .* sbet2(l) + dbet1 .* dbet2));

    if e2 ~= 0
      c2 = a^2 * (1 + (1 - e2) * eatanhe(1, e2) / e2) / 2;
    else
      c2 = a^2;
    end
    S12 = S12 + c2 * gam12;
    S12 = reshape(S12, S);
  end
end

function A1m1 = A1m1f(epsi)
%A1M1F  Evaluate A_1 - 1
%
%   A1m1 = A1M1F(epsi) evaluates A_1 - 1 using Eq. (17).  epsi and A1m1 are
%   K x 1 arrays.

  eps2 = epsi.^2;
  t = eps2.*(eps2.*(eps2+4)+64)/256;
  A1m1 = (t + epsi) ./ (1 - epsi);
end

function d = AngDiff(x, y)
%ANGDIFF  Compute angle difference accurately
%
%   d = ANGDIFF(x, y) computes y - x, reduces the result to (-180,180] and
%   rounds the result.  x and y must be in [-180,180].  x and y can be any
%   compatible shapes.

  [d, t] = sumx(-x, y);
  c = (d - 180) + t > 0;
  d(c) = (d(c) - 360) + t(c);
  c = (d + 180) + t <= 0;
  d(c) = (d(c) + 360) + t(c);
end

function x = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range [-180, 180)
%
%   x = ANGNORMALIZE(x) reduces angles in [-540, 540) to the range
%   [-180, 180).  x can be any shape.

  x(x >= 180) = x(x >= 180) - 360;
  x(x < -180) = x(x < -180) + 360;
end

function y = AngRound(x)
%ANGROUND  Round tiny values so that tiny values become zero.
%
%   y = ANGROUND(x) rounds x by adding and subtracting 1/16 to it if it is
%   small.  x can be any shape.

  z = 1/16;
  y = abs(x);
  y(y < z) = z - (z - y(y < z));
  y(x < 0) = -y(x < 0);
end

function C4 = C4f(epsi, C4x)
%C4F  Evaluate C_4
%
%   C4 = C4F(epsi, C4x) evaluates C_{4,l} in the expansion for the area
%   (Eq. (65) expressed in terms of n and epsi) using the coefficient
%   vector C4x.  epsi is a K x 1 array.  C4x is a 1 x 21 array.  C4 is a
%   K x 6 array.

  nC4 = 6;
  nC4x = size(C4x, 2);
  j = nC4x;
  C4 = zeros(length(epsi), nC4);
  for k = nC4 : -1 : 1
    t = C4(:, k);
    for i = nC4 - k : -1 : 0
      t = epsi .* t + C4x(j);
      j = j - 1;
    end
    C4(:, k) = t;
  end
  mult = ones(length(epsi), 1);
  for k = 2 : nC4
    mult = mult .* epsi;
    C4(:, k) = C4(:, k) .* mult;
  end
end

function ellipsoid = defaultellipsoid
%DEFAULTELLIPSOID  Return the WGS84 ellipsoid
%
%   ellipsoid = DEFAULTELLIPSOID
%
%   returns a vector of the equatorial radius and eccentricity for the
%   WGS84 ellipsoid.  use ecc2flat and flat2ecc to convert between
%   the eccentricity and the flattening.
%
%   See also ECC2FLAT, FLAT2ECC.

  a = 6378137;
  f = 1/298.257223563;
  e = flat2ecc(f);
  ellipsoid = [a, e];
end

function G4x = G4coeff(n)
%G4COEFF  Evaluate coefficients for C_4 for great ellipse
%
%   G4x = G4COEFF(n) evaluates the coefficients of epsilon^l in expansion
%   of the greate ellipse area (expressed in terms of n and epsi).  n is a
%   scalar.  G4x is a 1 x 21 array.

  nG4 = 6;
  nG4x = (nG4 * (nG4 + 1)) / 2;
  G4x = zeros(1, nG4x);
  G4x(0+1) = (n*(n*(n*(n*(200*n+416)+1144)+6864)+21021)+15015)/90090;
  G4x(1+1) = (n*(n*((-117944*n-110552)*n-84227)-41184)-9009)/120120;
  G4x(2+1) = (n*(n*(6417449*n+3013374)+1012583)+172458)/720720;
  G4x(3+1) = ((-135037988*n-32774196)*n-4232371)/5765760;
  G4x(4+1) = (138833443*n+13938873)/5765760;
  G4x(5+1) = -13200233/1537536;
  G4x(6+1) = (n*(n*(n*(117944*n+110552)+84227)+41184)+9009)/1081080;
  G4x(7+1) = (n*((-5975241*n-2676466)*n-847847)-136422)/4324320;
  G4x(8+1) = (n*(71379996*n+16424252)+1987557)/17297280;
  G4x(9+1) = (-39452953*n-3753828)/8648640;
  G4x(10+1) = 2625577/1537536;
  G4x(11+1) = (n*(n*(203633*n+80106)+20735)+2574)/1441440;
  G4x(12+1) = ((-3634676*n-741988)*n-76219)/5765760;
  G4x(13+1) = (2443153*n+208182)/2882880;
  G4x(14+1) = -5512967/15375360;
  G4x(15+1) = (n*(48020*n+8372)+715)/1153152;
  G4x(16+1) = (-71477*n-5317)/768768;
  G4x(17+1) = 22397/439296;
  G4x(18+1) = (1407*n+91)/329472;
  G4x(19+1) = -5453/1317888;
  G4x(20+1) = 21/146432;
end

function [x, y] = norm2(x, y)
%NORM2  Normalize x and y
%
%   [x, y] = NORM2(x, y) normalize x and y so that x^2 + y^2 = 1.  x and y
%   can be any shape.

  r = hypot(x, y);
  x = x ./ r;
  y = y ./ r;
end

function y = SinCosSeries(sinp, sinx, cosx, c)
%SINSCOSERIES  Evaluate a sine or cosine series using Clenshaw summation
%
%   y = SINCOSSERIES(sinp, sinx, cosx, c) evaluate
%     y = sum(c[i] * sin( 2*i    * x), i, 1, n), if  sinp
%     y = sum(c[i] * cos((2*i-1) * x), i, 1, n), if ~sinp
%
%   where n is the size of c.  x is given via its sine and cosine in sinx
%   and cosx.  sinp is a scalar.  sinx, cosx, and y are K x 1 arrays.  c is
%   a K x N array.

  if isempty(sinx), y = []; return, end
  n = size(c, 2);
  ar = 2 * (cosx - sinx) .* (cosx + sinx);
  y1 = zeros(length(sinx), 1);
  if mod(n, 2)
    y0 = c(:, n);
    n = n - 1;
  else
    y0 = y1;
  end

  for k = n : -2 : 1
    y1 = ar .* y0 - y1 + c(:, k);
    y0 = ar .* y1 - y0 + c(:, k-1);
  end
  if sinp
    y = 2 * sinx .* cosx .* y0;
  else
    y = cosx .* (y0 - y1);
  end
end

function [s, t] = sumx(u, v)
%SUM   Error free sum
%
%   [s, t] = SUMX(u, v) returns the rounded sum u + v in s and the error in
%   t, such that s + t = u + v, exactly.  u and v can be any compatible
%   shapes.

  s = u + v;
  up = s - v;
  vpp = s - up;
  up = up - u;
  vpp = vpp -  v;
  t = -(up + vpp);
end
function C1 = C1f(epsi)
%C1F  Evaluate C_{1,k}
%
%   C1 = C1F(epsi) evaluates C_{1,l} using Eq. (18).  epsi is a K x 1
%   array and C1 is a K x 6 array.

  nC1 = 6;
  C1 = zeros(length(epsi), nC1);
  eps2 = epsi.^2;
  d = epsi;
  C1(:,1) = d.*((6-eps2).*eps2-16)/32;
  d = d.*epsi;
  C1(:,2) = d.*((64-9*eps2).*eps2-128)/2048;
  d = d.*epsi;
  C1(:,3) = d.*(9*eps2-16)/768;
  d = d.*epsi;
  C1(:,4) = d.*(3*eps2-5)/512;
  d = d.*epsi;
  C1(:,5) = -7*d/1280;
  d = d.*epsi;
  C1(:,6) = -7*d/2048;
end
function y = eatanhe(x, e2)
%EATANHE   e*atanh(e*x)
%
%   EATANHE(x, e2) returns e*atanh(e*x) where e = sqrt(e2)
%   e2 is a scalar; x can be any shape.

  e = sqrt(abs(e2));
  if (e2 >= 0)
    y = e * atanh(e * x);
  else
    y = -e * atan(e * x);
  end
end