#include "vectorlib.h"

int main () {
  
  Vector a(3);
  for (int i=0; i<3; i++) {a(i) = i*i-1;};
  
  TimeVector A (3,4);
  for (int i=0; i<A.get_len(); i++) {A(i) = a+i;};
  
  a.print("a");
  A(2).print("A(2)");
  
  Vector b = a;
  b /= A(2);
  
  b.print("b");
  
  return 0;
};


