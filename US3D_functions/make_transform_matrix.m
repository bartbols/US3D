function T = make_transform_matrix(x,y,z,a,e,r,direction)

if nargin < 7
    direction = 'forward';
end

% Converted to Matlab from translate.cpp. This creates
% /* Pre-calculate radian angles */
deg2rad = pi/180;
ra = a*deg2rad;
rb = e*deg2rad;
rg = r*deg2rad;

% /* Pre-calculate sines and cosines */
cos_a = cos(ra);
cos_b = cos(rb);
cos_g = cos(rg);
sin_a = sin(ra);
sin_b = sin(rb);
sin_g = sin(rg);
%
m = zeros(16,1);
switch direction
        
    case 'forward'        
        m(1) = cos_a * cos_b;
        m(2) = sin_a * cos_b;
        m(3) = - sin_b;
        m(4) = 0.0;
        m(5) = cos_a * sin_b * sin_g - sin_a * cos_g;
        m(6) = sin_a * sin_b * sin_g + cos_a * cos_g;
        m(7) = cos_b * sin_g;
        m(8) = 0.0;
        m(9) = cos_a * sin_b * cos_g + sin_a * sin_g;
        m(10) = sin_a * sin_b * cos_g - cos_a * sin_g;
        m(11) = cos_b * cos_g;
        m(12) = 0.0;
        m(13) = x;
        m(14) = y;
        m(15) = z;
        m(16) = 1.0;
    case 'inverse'
        m(1) = cos_a * cos_b;
        m(2) = cos_a * sin_b * sin_g - sin_a * cos_g;
        m(3) = cos_a * sin_b * cos_g + sin_a * sin_g;
        m(4) = 0.0;
        m(5) = sin_a * cos_b;
        m(6) = sin_a * sin_b * sin_g + cos_a * cos_g;
        m(7) = sin_a * sin_b * cos_g - cos_a * sin_g;
        m(8) = 0.0;
        m(9) = - sin_b;
        m(10) = cos_b * sin_g;
        m(11) = cos_b * cos_g;
        m(12) = 0.0;
        m(13) = -x * m(1) -y * m(5) -z * m(9);
        m(14) = -x * m(2) -y * m(6) -z * m(10);
        m(15) = -x * m(3) -y * m(7) -z * m(11);
        m(16) = 1.0;
end

% Reshape into a 4x4 matrix.
T = reshape(m,4,4);