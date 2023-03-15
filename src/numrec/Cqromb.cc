//
// Complex version: can integrage complex# = function_of(real#)
//
// Romberg integration from Numerical Recipes
// (for the cautious at heart)
//
// NOTE: this has been modified to use a trpzd function without
//       static varriables.
//
#include <math.h>
#include "Globals.h"
#include "NRUtil.h"
#include "IEKG.h"

Real Cqromb(Real func(Real,Real*), Real a, Real b, Real eps, Real *args)
{
        int JMAX=100, JMAXP=JMAX+1, K=5;
        Complex ss,dss;
        Complex *s, *s_t;
        Real *h, *h_t;
        s = Complexvector(0,JMAX-1);
        h = Realvector(0,JMAXP-1);
        s_t = Complexvector(0,K-1);
        h_t = Realvector(0,K-1);
        int i,j;
                                                                                                                        
        h[0]=1.0;
        for (j=1;j<=JMAX;j++) {
                if(j==1){
                  Ctrpzd(func,a,b,j,0.0,&s[0],args);
                } else {
                  Ctrpzd(func,a,b,j,s[j-2],&s[j-1],args);
                }
                if (j >= K) {
                        for (i=0;i<K;i++) {
                                h_t[i]=h[j-K+i];
                                s_t[i]=s[j-K+i];
                        }
                        Cpolint(h_t,s_t,0.0,ss,dss,K);
                        //if (fabs(dss) <= eps*fabs(ss)) return ss;
                        //modified from above to account for complex 
                        //numbers and 0 result
                        if (fabs(dss.real()) <= eps*fabs(ss.real) &&
                            fabs(dss.imag()) <= eps*fabs(ss.imag) || 
                            fabs(ss.real()) < eps && fabs(ss.imag()) < eps) return ss;
                }
                h[j]=0.25*h[j-1];
        }
        Die("Too many steps in routine qromb");
        return 0.0;
}


