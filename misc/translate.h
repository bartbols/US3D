#ifndef translate_h
#define translate_h

#include <vector>
#define _USE_MATH_DEFINES 
#include <math.h>

// Matrix class for 4 by 4 transformation matrices
class matrix16d {
public:

  // Initialisation routines
  // default
  matrix16d() {}
  // copy from previous matrix
  matrix16d(double *m1);
  // create identity or zero matrix
  matrix16d(int i) {if (i==1) identity(); else zero(); }
  // create forward or inverse translation matrix
  matrix16d(double x, double y, double z, int i=1) { 
    if (i==-1) inverse(x, y, z); 
    else forward(x, y, z); }
  // create forward or inverse translation and rotation matrix using Euler angles
  matrix16d(double x, double y, double z, double a, double e, double r, int i=1) { 
    if (i==-1) inverse(x, y, z, a, e, r); 
    else forward(x, y, z, a, e, r); }

  // Functions which affect the current matrix
  // Set matrix to identity
  matrix16d& identity();
  // Set matrix to zero
  matrix16d& zero();
  // Set matrix to forward translation or rotation
  matrix16d& forward(double x, double y, double z);
  matrix16d& forward(double x, double y, double z, double a, double e, double r);
  // Set matrix to inverse translation or rotation
  matrix16d& inverse(double x, double y, double z);
  matrix16d& inverse(double x, double y, double z, double a, double e, double r);

  // Functions which return another matrix, leaving the current one unaffected
  // Generate the inverse of a transformation matrix
  matrix16d invert();
  // Extract translations and Euler rotations from matrix
  void decompose(double *x, double *y, double *z, double *a, double *e, double *r);

  // Matrix operations
  // * - multiply matrices
  matrix16d operator*(const matrix16d& m2);
  // *= - multiply and update matrix
  matrix16d& operator*=(const matrix16d& m2) { (*this)=(*this)*m2; return *this; }

  // Matrix component access
  // [i] - return element i of matrix 
  double& operator[](int i) { return m[i]; }
  // Expose a pointer to a set of doubles - can be used to initialise opengl matrices
  double m[16];

private:
};

// Vertex class for 3 component positon vectors
class vertex3d {
public:

  // Initialisation routines
  // default initialise to zero
  vertex3d() {x=0.0; y=0.0; z=0.0;}
  // initialise to given component values
  vertex3d(double a, double b, double c=0.0) {x=a; y=b; z=c;}
  // initialise using integers
  vertex3d(int a, int b, int c=0) {x=double(a); y=double(b); z=double(c);}

  // Vertex operators
  // == - test for equality
  bool operator==(const vertex3d& v) { 
    if ((x==v.x)&&(y==v.y)&&(z==v.z)) return true; else return false; }
  // != - test for in-equality
  bool operator!=(const vertex3d& v) { 
    if ((x!=v.x)||(y!=v.y)||(z!=v.z)) return true; else return false; }
  // + - vector addition
  vertex3d operator+(const vertex3d& v) { return vertex3d(x+v.x, y+v.y, z+v.z); }
  // - - vector subtraction
  vertex3d operator-(const vertex3d& v) { return vertex3d(x-v.x, y-v.y, z-v.z); }
  // += - vector addition and update
  vertex3d& operator+=(const vertex3d& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
  // -= - vector subtraction and update
  vertex3d& operator-=(const vertex3d& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
  // ^ - vector cross product, returning a vector
  vertex3d operator^(const vertex3d& v) { 
    return vertex3d(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
  // * - vector dot product, returning a double
  double operator*(const vertex3d& v) { return (x*v.x + y*v.y +z*v.z); }
  // * - vector multipled by a scalar (order not important)
  friend vertex3d operator*(const vertex3d& v, const double& a) { 
    return vertex3d(v.x*a, v.y*a, v.z*a); }
  friend vertex3d operator*(const double& a, const vertex3d& v) { 
    return vertex3d(v.x*a, v.y*a, v.z*a); }
  // *= - vector multiplied by scalar, updating the vector
  vertex3d& operator*=(const double& a) { x*=a; y*=a; z*=a; return *this; }
  // / - vector divided by scalar
  vertex3d operator/(const double& a) { return vertex3d(x/a, y/a, z/a); }
  // /= - vector divided by scalar, updating the scalar
  vertex3d& operator/=(const double& a) { x/=a; y/=a; z/=a; return *this; }

  // Functions not affecting the current vector
  // vector squared magnitude
  double abs2() { return (x*x + y*y + z*z); }
  // vector magnitude
  double abs() { return sqrt(this->abs2()); }
  // produce normalised version of the vector
  vertex3d norm() { double s(this->abs()); if (s==0) return *this; else return vertex3d(x/s, y/s, z/s); }
  // scale the vector by different values for each component
  vertex3d scale(double xs, double ys) { return vertex3d(x*xs, y*ys, z); }
  vertex3d scale(double xs, double ys, double zs) { return vertex3d(x*xs, y*ys, z*zs); }

  // * - matrix multiplied into a vector, returning a vector 
  friend vertex3d operator*(const matrix16d& m1, const vertex3d& v1);

  // vector component access
  double x, y, z;
private:
};

#endif /* translate_h */



















