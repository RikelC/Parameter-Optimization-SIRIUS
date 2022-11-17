#ifndef resolution_strip_h
#define resolution_strip_1
class resolution_strip{
  public:
  resolution_strip(unsigned short b, unsigned short c, unsigned short k1, unsigned short m1, double r){
    board = b;
    channel = c;
    k = k1;
    m = m1;
    FWHM = r;
  };
  unsigned short board;
  unsigned short channel;
  unsigned short k;
  unsigned short m;
  double FWHM;
};
#endif
