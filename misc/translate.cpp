#include "translate.h"

// Initialisation by copying from an existing matrix
matrix16d::matrix16d(double *m1) {
  m[0] = m1[0];
  m[1] = m1[1];
  m[2] = m1[2];
  m[3] = m1[3];
  m[4] = m1[4];
  m[5] = m1[5];
  m[6] = m1[6];
  m[7] = m1[7];
  m[8] = m1[8];
  m[9] = m1[9];
  m[10] = m1[10];
  m[11] = m1[11];
  m[12] = m1[12];
  m[13] = m1[13];
  m[14] = m1[14];
  m[15] = m1[15];
}

// Set entire matrix to zero
matrix16d& matrix16d::zero( )
{
  m[0] = m[4] = m[8] = m[12] = 0.0;
  m[1] = m[2] = m[3] = 0.0;
  m[5] = m[6] = m[7] = 0.0;
  m[9] = m[10] = m[11] = 0.0;
  m[13] = m[14] = m[15] = 0.0;
  return *this;
}

// Set to identity matrix
matrix16d& matrix16d::identity( )
{
  m[0] = m[5] = m[10] = m[15] = 1.0;
  m[1] = m[2] = m[3] = 0.0;
  m[4] = m[6] = m[7] = 0.0;
  m[8] = m[9] = m[11] = 0.0;
  m[12] = m[13] = m[14] = 0.0;
  return *this;
}

// Create forward translation only (no rotation)
matrix16d& matrix16d::forward(double x, double y, double z)
{
  m[0] = m[5] = m[10] = 1.0;
  m[1] = m[2] = m[3] = 0.0;
  m[4] = m[6] = m[7] = 0.0;
  m[8] = m[9] = m[11] = 0.0;
  m[12] = x;
  m[13] = y;
  m[14] = z;
  m[15] = 1.0;
  return *this;
}

// Create forward translation and rotation from Euler angles
// Note that this is actually Tait-Bryan in ZYX notation, with:
// alpha = azimuth = yaw = a
// beta = elevation = pitch = e
// gamma = roll = r
matrix16d& matrix16d::forward( double x, double y, double z, double a, double e, double r )
{
  double sin_a, sin_b, sin_g, cos_a, cos_b, cos_g;
  double ra, rb, rg;

  /* Pre-calculate radian angles */
  ra = a*M_PI/(double)180;
  rb = e*M_PI/(double)180;
  rg = r*M_PI/(double)180;

  /* Pre-calculate sines and cosines */
  cos_a = cos(ra);
  cos_b = cos(rb);
  cos_g = cos(rg);
  sin_a = sin(ra);
  sin_b = sin(rb);
  sin_g = sin(rg);

  /* Create the correct matrix coefficients */
  m[0] = cos_a * cos_b;
  m[1] = sin_a * cos_b;
  m[2] = - sin_b;
  m[3] = 0.0;
  m[4] = cos_a * sin_b * sin_g - sin_a * cos_g;
  m[5] = sin_a * sin_b * sin_g + cos_a * cos_g;
  m[6] = cos_b * sin_g;
  m[7] = 0.0;
  m[8] = cos_a * sin_b * cos_g + sin_a * sin_g;
  m[9] = sin_a * sin_b * cos_g - cos_a * sin_g;
  m[10] = cos_b * cos_g;
  m[11] = 0.0;
  m[12] = x;
  m[13] = y;
  m[14] = z;
  m[15] = 1.0;

  return *this;
}

// Create inverse translation only
matrix16d& matrix16d::inverse(double x, double y, double z)
{
  m[0] = m[5] = m[10] = 1.0;
  m[1] = m[2] = m[3] = 0.0;
  m[4] = m[6] = m[7] = 0.0;
  m[8] = m[9] = m[11] = 0.0;
  m[12] = -x;
  m[13] = -y;
  m[14] = -z;
  m[15] = 1.0;
  return *this;
}

// Create inverse translation and rotation from Euler angles 
// Note that this is actually Tait-Bryan in ZYX notation, with:
// alpha = azimuth = yaw = a
// beta = elevation = pitch = e
// gamma = roll = r
matrix16d& matrix16d::inverse( double x, double y, double z, double a, double e, double r )
{
  double sin_a, sin_b, sin_g, cos_a, cos_b, cos_g;
  double ra, rb, rg;

  /* Pre-calculate radian angles */
  ra = a*M_PI/(double)180;
  rb = e*M_PI/(double)180;
  rg = r*M_PI/(double)180;

  /* Pre-calculate sines and cosines */
  cos_a = cos(ra);
  cos_b = cos(rb);
  cos_g = cos(rg);
  sin_a = sin(ra);
  sin_b = sin(rb);
  sin_g = sin(rg);

  /* Create the correct matrix coefficients */
  m[0] = cos_a * cos_b;
  m[1] = cos_a * sin_b * sin_g - sin_a * cos_g;
  m[2] = cos_a * sin_b * cos_g + sin_a * sin_g;
  m[3] = 0.0;
  m[4] = sin_a * cos_b;
  m[5] = sin_a * sin_b * sin_g + cos_a * cos_g;
  m[6] = sin_a * sin_b * cos_g - cos_a * sin_g;
  m[7] = 0.0;
  m[8] = - sin_b;
  m[9] = cos_b * sin_g;
  m[10] = cos_b * cos_g;
  m[11] = 0.0;
  m[12] = -x * m[0] -y * m[4] -z * m[8];
  m[13] = -x * m[1] -y * m[5] -z * m[9];
  m[14] = -x * m[2] -y * m[6] -z * m[10];
  m[15] = 1.0;

  return *this;
}

// Specific inversion only suitable for 4 by 4 transformation matrices
matrix16d matrix16d::invert( )
{
  double zero_three, one_three, two_three;
  matrix16d mout(1);
  zero_three = -m[12]*m[0] - m[13]*m[1] - m[14]*m[2];
  one_three  = -m[12]*m[4] - m[13]*m[5] - m[14]*m[6];
  two_three  = -m[12]*m[8] - m[13]*m[9] - m[14]*m[10];
  mout.m[1] = m[4];
  mout.m[4] = m[1];
  mout.m[2] = m[8];
  mout.m[8] = m[2];
  mout.m[6] = m[9];
  mout.m[9] = m[6];
  mout.m[12] = zero_three;
  mout.m[13] = one_three;
  mout.m[14] = two_three;
  mout.m[0] = m[0];
  mout.m[5] = m[5];
  mout.m[10] = m[10];
  return mout;
}

// Calculate translation and Euler angles from 4 by 4 matrix
// Note that this is actually Tait-Bryan in ZYX notation, with:
// alpha = azimuth = yaw = a
// beta = elevation = pitch = e
// gamma = roll = r
void matrix16d::decompose( double *x, double *y, double *z, double *a, double *e, double *r )
{
  double tmp;

  // Translations are easy
  *x = m[12];
  *y = m[13];
  *z = m[14];

  // Catch degenerate elevation cases
  if (m[2] < -0.99999999) {
    *e = 90.0;
    *a = 0.0;
    *r = acos(m[8]);
    if ( (sin(*r)>0.0) ^ (m[4]>0.0) ) *r = -*r;
    *r *= 180.0/M_PI;
    return;
  }
  if (m[2] > 0.99999999) {
    *e = -90.0;
    *a = 0.0;
    *r = acos(m[5]);
    if ( (sin(*r)<0.0) ^ (m[4]>0.0) ) *r = -*r;
    *r *= 180.0/M_PI;
    return;
  }

  // Non-degenerate elevation - between -90 and +90
  *e = asin( -m[2] );

  // Now work out azimuth - between -180 and +180
  tmp = m[0]/cos(*e); // the denominator will not be zero
  if ( tmp <= -1.0 ) *a = M_PI;
  else if ( tmp >= 1.0 ) *a = 0.0;
  else *a = acos( tmp );
  if ( ((sin(*a) * cos(*e))>0.0) ^ ((m[1])>0.0) ) *a = -*a;

  // Now work out roll - between -180 and +180
  tmp = m[10]/cos(*e); // the denominator will not be zero
  if ( tmp <= -1.0 ) *r = M_PI;
  else if ( tmp >= 1.0 ) *r = 0.0;
  else *r = acos( tmp );
  if ( ((sin(*r) * cos(*e))>0.0) ^ ((m[6])>0.0) ) *r = -*r;

  // The output angles are in degrees
  *e *= 180.0/M_PI;
  *a *= 180.0/M_PI;
  *r *= 180.0/M_PI;
}

// Specific multiplication only suitable for transformation matrices
matrix16d matrix16d::operator*(const matrix16d& m2)
     /* multiply calibration matrix into position matrix */
     /* first = first * second                           */
{
  matrix16d mout(1);

  /* Note that we know of the existence of some zeros and ones */
  mout.m[0] = m[0] * m2.m[0] + m[4] * m2.m[1] + m[8] * m2.m[2];
  mout.m[1] = m[1] * m2.m[0] + m[5] * m2.m[1] + m[9] * m2.m[2];
  mout.m[2] = m[2] * m2.m[0] + m[6] * m2.m[1] + m[10] * m2.m[2];
  mout.m[4] = m[0] * m2.m[4] + m[4] * m2.m[5] + m[8] * m2.m[6];
  mout.m[5] = m[1] * m2.m[4] + m[5] * m2.m[5] + m[9] * m2.m[6];
  mout.m[6] = m[2] * m2.m[4] + m[6] * m2.m[5] + m[10] * m2.m[6];
  mout.m[8] = m[0] * m2.m[8] + m[4] * m2.m[9] + m[8] * m2.m[10];
  mout.m[9] = m[1] * m2.m[8] + m[5] * m2.m[9] + m[9] * m2.m[10];
  mout.m[10] = m[2] * m2.m[8] + m[6] * m2.m[9] + m[10] * m2.m[10];
  mout.m[12] = m[0] * m2.m[12] + m[4] * m2.m[13] + m[8] * m2.m[14] + m[12];
  mout.m[13] = m[1] * m2.m[12] + m[5] * m2.m[13] + m[9] * m2.m[14] + m[13];
  mout.m[14] = m[2] * m2.m[12] + m[6] * m2.m[13] + m[10] * m2.m[14] + m[14];
  return mout;
}

// Multiple vector by a transformation matrix
vertex3d operator*(const matrix16d& m1, const vertex3d& v1)
{
  vertex3d v;
  if ( v1.z == 0.0 ) {
    v.x = v1.x * m1.m[0] + v1.y * m1.m[4] + m1.m[12];
    v.y = v1.x * m1.m[1] + v1.y * m1.m[5] + m1.m[13];
    v.z = v1.x * m1.m[2] + v1.y * m1.m[6] + m1.m[14];
  } else {
    v.x = v1.x * m1.m[0] + v1.y * m1.m[4] + v1.z * m1.m[8]  + m1.m[12];
    v.y = v1.x * m1.m[1] + v1.y * m1.m[5] + v1.z * m1.m[9]  + m1.m[13];
    v.z = v1.x * m1.m[2] + v1.y * m1.m[6] + v1.z * m1.m[10] + m1.m[14];
  }
  return v;
}

