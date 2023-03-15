//
// Complex version of trpzd
// 
// NOTE: This has been modified so as not to use static varriables.
//

#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

void Ctrpzd(Complex func(Real,Real*), Real a, Real b, int n, Complex last, Complex *next, Real *args)
{
        Real x,tnm,sum,del;
        int it,j;

        if (n == 1) {
                *next=0.5*(b-a)*(func(a,args)+func(b,args));
                return;
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x,args);
                *next=0.5*(last+(b-a)*sum/tnm);
                return;
        }
}
