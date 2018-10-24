#ifndef VECTORLIB_H
#define VECTORLIB_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

class Vector;
class TimeVector;

class Vector {
  private:
    int     dim;
    double* data;
  public:
    // constructors and destructor
    Vector     ();
    void init  (int n, double c=0);
    Vector     (int n, double c=0);
    Vector     (const Vector& copy);
    ~Vector    ();
    // misc
    int     get_dim    () const;
    double* get_vector () const;
    
    void print      (string name="vector", int precision=3) const;
    void print_norm (string name="vector", int precision=3) const;
    // overloaded operators
    double  operator () (int i) const;
    double& operator () (int i);
    Vector& operator =  (const Vector& V);
    Vector& operator += (const Vector& V);
    Vector& operator -= (const Vector& V);
    Vector& operator *= (const Vector& V);
    Vector& operator /= (const Vector& V);
    Vector& operator += (double c);
    Vector& operator -= (double c);
    Vector& operator *= (double c);
    Vector& operator /= (double c);
    // additional operations
    double norm   ()                const;
    double norm_2 ()                const;
    double dot    (const Vector& V) const;
};

class TimeVector {
  private:
    int     dim;
    int     N;
    Vector* timedata;
  public:
    // constructors and destructors
    TimeVector  ();
    void init   (int d, int n);
    TimeVector  (int d, int n);
    ~TimeVector ();
    // misc
    int get_dim () const;
    int get_len () const;
    // overloaded operators
    Vector  operator () (int i) const;
    Vector& operator () (int i);
    // additional operations
    double norm (int i) const;
};

// vector operations
Vector operator - (const Vector& V);
Vector operator + (const Vector& V1, const Vector& V2);
Vector operator - (const Vector& V1, const Vector& V2);
Vector operator * (const Vector& V1, const Vector& V2);
Vector operator / (const Vector& V1, const Vector& V2);
Vector operator % (const Vector& V1, const Vector& V2);  // vector cross-product, assumes dim=3
// scalar operations
Vector operator + (const Vector& V, double c);
Vector operator - (const Vector& V, double c);
Vector operator * (const Vector& V, double c);
Vector operator / (const Vector& V, double c);
Vector operator + (double c, const Vector& V);
Vector operator - (double c, const Vector& V);
Vector operator * (double c, const Vector& V);
Vector operator / (double c, const Vector& V);


#endif
