#include "vectorlib.h"

/*
-------------------------------------
            Vector Class
--------------------------------------
*/

Vector::Vector () {
  dim  = 0;
};

void Vector::init (int n, double c) {
  dim  = n;
  data = new double [n];
  for (int i=0; i<n; i++) {data[i] = c;};
};

Vector::Vector (int n, double c) {
  dim  = n;
  data = new double [n];
  for (int i=0; i<n; i++) {data[i] = c;};
};

Vector::Vector (const Vector& copy) {
  dim  = copy.get_dim();
  data = new double [dim];
  for (int i=0; i<dim; i++) {data[i] = copy.data[i];};
};

Vector::~Vector () {
  delete[] data;
  dim  = 0;
  data = NULL;
};

int     Vector::get_dim    () const {return dim;};
double* Vector::get_vector () const {return data;};

double  Vector::operator () (int i) const     {return data[i];};
double& Vector::operator () (int i)           {return data[i];};

Vector& Vector::operator = (const Vector& V) {
  dim = V.get_dim();
  for (int i=0; i<dim; i++) {data[i] = V.data[i];};
  return *this;
};

Vector& Vector::operator += (const Vector& V) {
  if (dim == V.get_dim()) {for (int i=0; i<dim; i++) {data[i] += V.data[i];};};
  return *this;
};

Vector& Vector::operator -= (const Vector& V) {
  if (dim == V.get_dim()) {for (int i=0; i<dim; i++) {data[i] -= V.data[i];};};
  return *this;
};

Vector& Vector::operator *= (const Vector& V) {
  if (dim == V.get_dim()) {for (int i=0; i<dim; i++) {data[i] *= V.data[i];};};
  return *this;
};

Vector& Vector::operator /= (const Vector& V) {
  if (dim == V.get_dim()) {for (int i=0; i<dim; i++) {data[i] /= V.data[i];};};
  return *this;
};

void Vector::print (string name, int precision) const {
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << showpos << setprecision(precision) << scientific;
  cout << name    << " = " <<  endl;
  for (int i=0; i<dim; i++) {cout << "  " << data[i] << endl;};
  cout.flags(f);
};

void Vector::print_norm (string name, int precision) const {
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << showpos << setprecision(precision) << scientific;
  cout << "|"     << name  << "| = "         << norm() << endl;
  cout.flags(f);
};

double Vector::norm () const {
  double norm = 0.0;
  for (int i=0; i<dim; i++) {norm += data[i]*data[i];};
  return pow(norm,0.5);
};

/*
-------------------------------------
       Vector Class Operations
--------------------------------------
*/

Vector operator - (const Vector& V) {
  Vector negV(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {negV(i) = -V(i);};
  return negV;
};

Vector operator + (const Vector& V1, const Vector& V2) {
  if (V1.get_dim() == V2.get_dim()) {
    Vector V3(V1.get_dim());
    for (int i=0; i<V1.get_dim(); i++) {V3(i) = V1(i)+V2(i);};
    return V3;
  } else {return 0;};
};

Vector operator - (const Vector& V1, const Vector& V2) {
  if (V1.get_dim() == V2.get_dim()) {
    Vector V3(V1.get_dim());
    for (int i=0; i<V1.get_dim(); i++) {V3(i) = V1(i)-V2(i);};
    return V3;
  } else {return 0;};
};

Vector operator * (const Vector& V1, const Vector& V2) {
  if (V1.get_dim() == V2.get_dim()) {
    Vector V3(V1.get_dim());
    for (int i=0; i<V1.get_dim(); i++) {V3(i) = V1(i)*V2(i);};
    return V3;
  } else {return 0;};
};

Vector operator / (const Vector& V1, const Vector& V2) {
  if (V1.get_dim() == V2.get_dim()) {
    Vector V3(V1.get_dim());
    for (int i=0; i<V1.get_dim(); i++) {V3(i) = V1(i)/V2(i);};
    return V3;
  } else {return 0;};
};

Vector operator + (const Vector& V, double c) {
  Vector Vpc(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {Vpc(i) = V(i)+c;};
  return Vpc;
};

Vector operator + (double c, const Vector& V) {return V+c;};

Vector operator - (const Vector& V, double c) {
  Vector Vmc(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {Vmc(i) = V(i)-c;};
  return Vmc;
};

Vector operator - (double c, const Vector& V) {
  Vector cmV(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {cmV(i) = c-V(i);};
  return cmV;
};

Vector operator * (const Vector& V, double c) {
  Vector cxV(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {cxV(i) = c*V(i);};
  return cxV;
};

Vector operator * (double c, const Vector& V) {return V*c;};

Vector operator / (const Vector& V, double c) {
  Vector Vdc(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {Vdc(i) = V(i)/c;};
  return Vdc;
};

Vector operator / (double c, const Vector& V) {
  Vector cdV(V.get_dim());
  for (int i=0; i<V.get_dim(); i++) {cdV(i) = c/V(i);};
  return cdV;
};

/*
-------------------------------------
          TimeVector Class
--------------------------------------
*/

TimeVector::TimeVector () {
  dim      = 0;
  N        = 0;
};

void TimeVector::init (int d, int n) {
  dim      = d;
  N        = n;
  timedata = new Vector [n];
  for (int i=0; i<n; i++) {timedata[i].init(d);};
};

TimeVector::TimeVector (int d, int n) {
  dim      = d;
  N        = n;
  timedata = new Vector [n];
  for (int i=0; i<n; i++) {timedata[i].init(d);};
};

TimeVector::~TimeVector () {
  delete[] timedata;
  timedata = NULL;
  dim      = 0;
  N        = 0;
};

Vector  TimeVector::operator () (int i) const {return timedata[i];};
Vector& TimeVector::operator () (int i)       {return timedata[i];};

int TimeVector::get_dim () const {return dim;};
int TimeVector::get_len () const {return N;};

double TimeVector::norm (int i) const {return timedata[i].norm();};

