#include "luclass.h"

void first(hls::vector<hls::vector<double,SIZE>,SIZE> &l, hls::vector<hls::vector<double,SIZE>,SIZE> &u)
{
    Matrix matrix;
    matrix.randfill();
    matrix.LUPivot(l, u);
}
