#include <R.h>

void DescendMin(double *yvals, int *numin, int *istart,
                int *ilower, int *iupper) {

    int i;

    for (i = *istart; i > 0; i--)
        if (yvals[i-1] >= yvals[i])
            break;
    *ilower = i;

    for (i = *istart; i < *numin-1; i++)
        if (yvals[i+1] >= yvals[i])
            break;
    *iupper = i;
}

void FindEqualGreaterUnsorted(const double *in, const int *size, const double *target,
                      int *index) {
   int i;

   for (i=0; (i< *size-1 &&  in[i] < *target); i++) {}

   *index = i;
}
