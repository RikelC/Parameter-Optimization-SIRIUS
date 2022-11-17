#ifndef trap_parameter_h
#define trap_parameter_1
class trap_parameter{
  public:
  trap_parameter( unsigned short a, unsigned short b, double c){
    k = a; m = b; FWHM = c;
  };
  unsigned short k;
  unsigned short m;
  double FWHM;
};
#endif
