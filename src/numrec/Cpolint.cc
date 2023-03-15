#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

void polint(Real *xa, Real *ya, Real x, Real &y, Real &dy, int n)
{
        int i,m,ns=0;
        Real den,dif,dift,ho,hp,w;
        Real *c, *d;
        c = Realvector(0,n-1);
        d = Realvector(0,n-1);
        dif=fabs(x-xa[0]);
        for (i=0;i<n;i++) {
                if ((dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        y=ya[ns--];
        for (m=1;m<n;m++) {
                for (i=0;i<n-m;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+m]-x;
                        w=c[i+1]-d[i];
                        if ((den=ho-hp) == 0.0) Die("Error in routine polint");
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
        }
}
