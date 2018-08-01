/// @file       c74_lib_fft.h
///	@ingroup 	minlib
///	@copyright	Copyright 2018 The Min-Lib Authors. All rights reserved.
///	@license	Use of this source code is governed by the MIT License found in the License.md file.

#pragma once

#include "c74_min_api.h"


namespace c74 { namespace min { namespace lib {
    
    /// Defines a few support classes containing messy FFT math. Currently all
    /// the code is taken from the max/c74support-private/fft-support folder and
    /// looks a bit harsh
    
    /// TODO: Clean fftSupport, make it look like clean C++ code
    
    namespace fftSupport{
        
        /// Static trig functions and variables for the FFT
        
        class fftTrig{
        public:
            
        #ifdef REAL
        #else
        #define REAL double
        #endif
            
        #ifdef GOOD_TRIG
        #else
        #define FAST_TRIG
        #endif
            
        #if defined(GOOD_TRIG)
            
        #define FHT_SWAP(a,b,t) {(t)=(a);(a)=(b);(b)=(t);}
            
        #define TRIG_VARS					         \
            int t_lam=0; \
            REAL coswrk[20], sinwrk[20];	\
            coswrk[0] = 0.; \
            coswrk[1] = .70710678118654752440084436210484903928483593768846;\
            sinwrk[0] = 1.;\
            sinwrk[1] = .70710678118654752440084436210484903928483593768846;
            
        #define TRIG_INIT(k,c,s)					 \
        {								 \
            int i;							 \
            for (i=2 ; i<=k ; i++)					 \
            {coswrk[i]=costab[i];sinwrk[i]=sintab[i];}		 \
            t_lam = 0;					         \
            c = 1;							 \
            s = 0;							 \
        }
            
        #define TRIG_NEXT(k,c,s)					 \
        {								 \
            int i,j;	                                         \
            (t_lam)++;	  					 \
            for (i=0 ; !((1<<i)&t_lam) ; i++);			 \
            i = k-i;						 \
            s = sinwrk[i];						 \
            c = coswrk[i];						 \
            if (i>1)   						 \
            {	    						 \
                for (j=k-i+2 ; (1<<j)&t_lam ; j++);		 \
                j	       = k - j;					 \
                sinwrk[i] = halsec[i] * (sinwrk[i-1] + sinwrk[j]);  \
                coswrk[i] = halsec[i] * (coswrk[i-1] + coswrk[j]);  \
            }                                                    \
        }
            
        #define TRIG_RESET(k,c,s)
        #endif
            
        #if defined(FAST_TRIG)
            
        #define TRIG_VARS					 \
            REAL t_c,t_s;
            
        #define TRIG_INIT(k,c,s)				 \
        {				        		 \
            t_c  = costab[k];				         \
            t_s  = sintab[k];				         \
            c    = 1;				    		 \
            s    = 0;				    		 \
        }
            
        #define TRIG_NEXT(k,c,s)				 \
        {                                                    \
            REAL t = c;                                         \
            c   = t*t_c - s*t_s;				 \
            s   = t*t_s + s*t_c;				 \
        }
            
        #define TRIG_RESET(k,c,s)
        #endif
            
            
            REAL halsec[20]=
            {
                0,
                0,
                .54119610014619698439972320536638942006107206337801,
                .50979557910415916894193980398784391368261849190893,
                .50241928618815570551167011928012092247859337193963,
                .50060299823519630134550410676638239611758632599591,
                .50015063602065098821477101271097658495974913010340,
                .50003765191554772296778139077905492847503165398345,
                .50000941253588775676512870469186533538523133757983,
                .50000235310628608051401267171204408939326297376426,
                .50000058827484117879868526730916804925780637276181,
                .50000014706860214875463798283871198206179118093251,
                .50000003676714377807315864400643020315103490883972,
                .50000000919178552207366560348853455333939112569380,
                .50000000229794635411562887767906868558991922348920,
                .50000000057448658687873302235147272458812263401372
            };
            
            REAL costab[20]=
            {
                .00000000000000000000000000000000000000000000000000,
                .70710678118654752440084436210484903928483593768847,
                .92387953251128675612818318939678828682241662586364,
                .98078528040323044912618223613423903697393373089333,
                .99518472667219688624483695310947992157547486872985,
                .99879545620517239271477160475910069444320361470461,
                .99969881869620422011576564966617219685006108125772,
                .99992470183914454092164649119638322435060646880221,
                .99998117528260114265699043772856771617391725094433,
                .99999529380957617151158012570011989955298763362218,
                .99999882345170190992902571017152601904826792288976,
                .99999970586288221916022821773876567711626389934930,
                .99999992646571785114473148070738785694820115568892,
                .99999998161642929380834691540290971450507605124278,
                .99999999540410731289097193313960614895889430318945,
                .99999999885102682756267330779455410840053741619428
            };
            
            REAL sintab[20]=
            {
                1.0000000000000000000000000000000000000000000000000,
                .70710678118654752440084436210484903928483593768846,
                .38268343236508977172845998403039886676134456248561,
                .19509032201612826784828486847702224092769161775195,
                .09801714032956060199419556388864184586113667316749,
                .04906767432741801425495497694268265831474536302574,
                .02454122852291228803173452945928292506546611923944,
                .01227153828571992607940826195100321214037231959176,
                .00613588464915447535964023459037258091705788631738,
                .00306795676296597627014536549091984251894461021344,
                .00153398018628476561230369715026407907995486457522,
                .00076699031874270452693856835794857664314091945205,
                .00038349518757139558907246168118138126339502603495,
                .00019174759731070330743990956198900093346887403385,
                .00009587379909597734587051721097647635118706561284,
                .00004793689960306688454900399049465887274686668768
            };
            
            REAL coswrk[20]=
            {
                .00000000000000000000000000000000000000000000000000,
                .70710678118654752440084436210484903928483593768847,
                .92387953251128675612818318939678828682241662586364,
                .98078528040323044912618223613423903697393373089333,
                .99518472667219688624483695310947992157547486872985,
                .99879545620517239271477160475910069444320361470461,
                .99969881869620422011576564966617219685006108125772,
                .99992470183914454092164649119638322435060646880221,
                .99998117528260114265699043772856771617391725094433,
                .99999529380957617151158012570011989955298763362218,
                .99999882345170190992902571017152601904826792288976,
                .99999970586288221916022821773876567711626389934930,
                .99999992646571785114473148070738785694820115568892,
                .99999998161642929380834691540290971450507605124278,
                .99999999540410731289097193313960614895889430318945,
                .99999999885102682756267330779455410840053741619428
            };
            
            REAL sinwrk[20]=
            {
                1.0000000000000000000000000000000000000000000000000,
                .70710678118654752440084436210484903928483593768846,
                .38268343236508977172845998403039886676134456248561,
                .19509032201612826784828486847702224092769161775195,
                .09801714032956060199419556388864184586113667316749,
                .04906767432741801425495497694268265831474536302574,
                .02454122852291228803173452945928292506546611923944,
                .01227153828571992607940826195100321214037231959176,
                .00613588464915447535964023459037258091705788631738,
                .00306795676296597627014536549091984251894461021344,
                .00153398018628476561230369715026407907995486457522,
                .00076699031874270452693856835794857664314091945205,
                .00038349518757139558907246168118138126339502603495,
                .00019174759731070330743990956198900093346887403385,
                .00009587379909597734587051721097647635118706561284,
                .00004793689960306688454900399049465887274686668768
            };
            
        };
        
        
        /// Class containing FFT and IFFT transforms
        
        class fftTransforms : public fftTrig {
        public:
        #define SQRT2   2*0.70710678118654752440084436210484
            void fht(float *fz, int n)
            {
                int i,k,k1,k2,k3,k4,kx;
                float /* REAL */ *fi,*fn,*gi;
                TRIG_VARS;
                
                for (k1=1,k2=0;k1<n;k1++) {
                    REAL a;
                    for (k=n>>1; (!((k2^=k)&k)); k>>=1);
                    if (k1>k2) {
                        a=fz[k1];
                        fz[k1]=fz[k2];
                        fz[k2]=a;
                    }
                }
                for ( k=0 ; (1<<k)<n ; k++ );
                k  &= 1;
                if (k==0) {
                    for (fi=fz,fn=fz+n;fi<fn;fi+=4) {
                        REAL f0,f1,f2,f3;
                        f1     = fi[0 ]-fi[1 ];
                        f0     = fi[0 ]+fi[1 ];
                        f3     = fi[2 ]-fi[3 ];
                        f2     = fi[2 ]+fi[3 ];
                        fi[2 ] = (f0-f2);
                        fi[0 ] = (f0+f2);
                        fi[3 ] = (f1-f3);
                        fi[1 ] = (f1+f3);
                    }
                } else {
                    for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi+=8,gi+=8) {
                        REAL s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
                        c1     = fi[0 ] - gi[0 ];
                        s1     = fi[0 ] + gi[0 ];
                        c2     = fi[2 ] - gi[2 ];
                        s2     = fi[2 ] + gi[2 ];
                        c3     = fi[4 ] - gi[4 ];
                        s3     = fi[4 ] + gi[4 ];
                        c4     = fi[6 ] - gi[6 ];
                        s4     = fi[6 ] + gi[6 ];
                        f1     = (s1 - s2);
                        f0     = (s1 + s2);
                        g1     = (c1 - c2);
                        g0     = (c1 + c2);
                        f3     = (s3 - s4);
                        f2     = (s3 + s4);
                        g3     = SQRT2*c4;
                        g2     = SQRT2*c3;
                        fi[4 ] = f0 - f2;
                        fi[0 ] = f0 + f2;
                        fi[6 ] = f1 - f3;
                        fi[2 ] = f1 + f3;
                        gi[4 ] = g0 - g2;
                        gi[0 ] = g0 + g2;
                        gi[6 ] = g1 - g3;
                        gi[2 ] = g1 + g3;
                    }
                }
                if (n<16)
                    return;
                
                do {
                    double /* REAL */ s1,c1;
                    k  += 2;
                    k1  = 1  << k;
                    k2  = k1 << 1;
                    k4  = k2 << 1;
                    k3  = k2 + k1;
                    kx  = k1 >> 1;
                    fi  = fz;
                    gi  = fi + kx;
                    fn  = fz + n;
                    do {
                        REAL g0,f0,f1,g1,f2,g2,f3,g3;
                        f1      = fi[0 ] - fi[k1];
                        f0      = fi[0 ] + fi[k1];
                        f3      = fi[k2] - fi[k3];
                        f2      = fi[k2] + fi[k3];
                        fi[k2]  = f0	  - f2;
                        fi[0 ]  = f0	  + f2;
                        fi[k3]  = f1	  - f3;
                        fi[k1]  = f1	  + f3;
                        g1      = gi[0 ] - gi[k1];
                        g0      = gi[0 ] + gi[k1];
                        g3      = SQRT2  * gi[k3];
                        g2      = SQRT2  * gi[k2];
                        gi[k2]  = g0	  - g2;
                        gi[0 ]  = g0	  + g2;
                        gi[k3]  = g1	  - g3;
                        gi[k1]  = g1	  + g3;
                        gi     += k4;
                        fi     += k4;
                    } while (fi<fn);
                    TRIG_INIT(k,c1,s1);
                    for (i=1;i<kx;i++) {
                        double /* REAL */ c2,s2;
                        TRIG_NEXT(k,c1,s1);
                        c2 = c1*c1 - s1*s1;
                        s2 = 2*(c1*s1);
                        fn = fz + n;
                        fi = fz +i;
                        gi = fz +k1-i;
                        do {
                            double /* REAL */ a,b,g0,f0,f1,g1,f2,g2,f3,g3;
                            b       = s2*fi[k1] - c2*gi[k1];
                            a       = c2*fi[k1] + s2*gi[k1];
                            f1      = fi[0 ]    - a;
                            f0      = fi[0 ]    + a;
                            g1      = gi[0 ]    - b;
                            g0      = gi[0 ]    + b;
                            b       = s2*fi[k3] - c2*gi[k3];
                            a       = c2*fi[k3] + s2*gi[k3];
                            f3      = fi[k2]    - a;
                            f2      = fi[k2]    + a;
                            g3      = gi[k2]    - b;
                            g2      = gi[k2]    + b;
                            b       = s1*f2     - c1*g3;
                            a       = c1*f2     + s1*g3;
                            fi[k2]  = f0        - a;
                            fi[0 ]  = f0        + a;
                            gi[k3]  = g1        - b;
                            gi[k1]  = g1        + b;
                            b       = c1*g2     - s1*f3;
                            a       = s1*g2     + c1*f3;
                            gi[k2]  = g0        - a;
                            gi[0 ]  = g0        + a;
                            fi[k3]  = f1        - b;
                            fi[k1]  = f1        + b;
                            gi     += k4;
                            fi     += k4;
                        } while (fi<fn);
                    }
                    TRIG_RESET(k,c1,s1);
                } while (k4<n);
            }
            
            void ifft(int n, float *real, float *imag)
            {
                double a,b,c,d;
                double q,r,s,t;
                int i,j,k;
                fht(real,n);
                fht(imag,n);
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i]; b = real[j];  q=a+b; r=a-b;
                    c = imag[i]; d = imag[j];  s=c+d; t=c-d;
                    imag[i] = (s+r)*0.5;  imag[j] = (s-r)*0.5;
                    real[i] = (q-t)*0.5;  real[j] = (q+t)*0.5;
                }
            }
            
            void realfft(int n, float *real)
            {
                double a,b;
                int i,j,k;
                
                fht(real,n);
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i];
                    b = real[j];
                    real[j] = (a-b)*0.5;
                    real[i] = (a+b)*0.5;
                }
            }
            
            void fft(int n, float *real, float *imag)
            {
                double a,b,c,d;
                double q,r,s,t;
                int i,j,k;
                
                #ifdef DBG
                debug_printf("fft real %.2f %.2f",
                             real[0],real[1]);
                
                debug_printf("fft imag %.2f %.2f",
                             imag[0],imag[1]);
                #endif
                
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i]; b = real[j];  q=a+b; r=a-b;
                    c = imag[i]; d = imag[j];  s=c+d; t=c-d;
                    real[i] = (q+t)*.5; real[j] = (q-t)*.5;
                    imag[i] = (s-r)*.5; imag[j] = (s+r)*.5;
                }
                fht(real,n);
                fht(imag,n);
                #ifdef DBG
                debug_printf("fft after real %.2f %.2f",
                             real[0],real[1]);
                
                debug_printf("fft after imag %.2f %.2f",
                             imag[0],imag[1]);
                #endif
            }
            
            void realifft(int n, float *real)
            {
                double a,b;
                int i,j,k;
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i];
                    b = real[j];
                    real[j] = (a-b);
                    real[i] = (a+b);
                }
                fht(real,n);
            }
            
            
            void fht64(double *fz, int n)
            {
                int i,k,k1,k2,k3,k4,kx;
                double /* REAL */ *fi,*fn,*gi;
                TRIG_VARS;
                
                for (k1=1,k2=0;k1<n;k1++) {
                    REAL a;
                    for (k=n>>1; (!((k2^=k)&k)); k>>=1);
                    if (k1>k2) {
                        a=fz[k1];
                        fz[k1]=fz[k2];
                        fz[k2]=a;
                    }
                }
                for ( k=0 ; (1<<k)<n ; k++ );
                k  &= 1;
                if (k==0) {
                    for (fi=fz,fn=fz+n;fi<fn;fi+=4) {
                        REAL f0,f1,f2,f3;
                        f1     = fi[0 ]-fi[1 ];
                        f0     = fi[0 ]+fi[1 ];
                        f3     = fi[2 ]-fi[3 ];
                        f2     = fi[2 ]+fi[3 ];
                        fi[2 ] = (f0-f2);
                        fi[0 ] = (f0+f2);
                        fi[3 ] = (f1-f3);
                        fi[1 ] = (f1+f3);
                    }
                } else {
                    for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi+=8,gi+=8) {
                        REAL s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
                        c1     = fi[0 ] - gi[0 ];
                        s1     = fi[0 ] + gi[0 ];
                        c2     = fi[2 ] - gi[2 ];
                        s2     = fi[2 ] + gi[2 ];
                        c3     = fi[4 ] - gi[4 ];
                        s3     = fi[4 ] + gi[4 ];
                        c4     = fi[6 ] - gi[6 ];
                        s4     = fi[6 ] + gi[6 ];
                        f1     = (s1 - s2);
                        f0     = (s1 + s2);
                        g1     = (c1 - c2);
                        g0     = (c1 + c2);
                        f3     = (s3 - s4);
                        f2     = (s3 + s4);
                        g3     = SQRT2*c4;
                        g2     = SQRT2*c3;
                        fi[4 ] = f0 - f2;
                        fi[0 ] = f0 + f2;
                        fi[6 ] = f1 - f3;
                        fi[2 ] = f1 + f3;
                        gi[4 ] = g0 - g2;
                        gi[0 ] = g0 + g2;
                        gi[6 ] = g1 - g3;
                        gi[2 ] = g1 + g3;
                    }
                }
                if (n<16)
                    return;
                
                do {
                    double /* REAL */ s1,c1;
                    k  += 2;
                    k1  = 1  << k;
                    k2  = k1 << 1;
                    k4  = k2 << 1;
                    k3  = k2 + k1;
                    kx  = k1 >> 1;
                    fi  = fz;
                    gi  = fi + kx;
                    fn  = fz + n;
                    do {
                        REAL g0,f0,f1,g1,f2,g2,f3,g3;
                        f1      = fi[0 ] - fi[k1];
                        f0      = fi[0 ] + fi[k1];
                        f3      = fi[k2] - fi[k3];
                        f2      = fi[k2] + fi[k3];
                        fi[k2]  = f0	  - f2;
                        fi[0 ]  = f0	  + f2;
                        fi[k3]  = f1	  - f3;
                        fi[k1]  = f1	  + f3;
                        g1      = gi[0 ] - gi[k1];
                        g0      = gi[0 ] + gi[k1];
                        g3      = SQRT2  * gi[k3];
                        g2      = SQRT2  * gi[k2];
                        gi[k2]  = g0	  - g2;
                        gi[0 ]  = g0	  + g2;
                        gi[k3]  = g1	  - g3;
                        gi[k1]  = g1	  + g3;
                        gi     += k4;
                        fi     += k4;
                    } while (fi<fn);
                    TRIG_INIT(k,c1,s1);
                    for (i=1;i<kx;i++) {
                        double /* REAL */ c2,s2;
                        TRIG_NEXT(k,c1,s1);
                        c2 = c1*c1 - s1*s1;
                        s2 = 2*(c1*s1);
                        fn = fz + n;
                        fi = fz +i;
                        gi = fz +k1-i;
                        do {
                            double /* REAL */ a,b,g0,f0,f1,g1,f2,g2,f3,g3;
                            b       = s2*fi[k1] - c2*gi[k1];
                            a       = c2*fi[k1] + s2*gi[k1];
                            f1      = fi[0 ]    - a;
                            f0      = fi[0 ]    + a;
                            g1      = gi[0 ]    - b;
                            g0      = gi[0 ]    + b;
                            b       = s2*fi[k3] - c2*gi[k3];
                            a       = c2*fi[k3] + s2*gi[k3];
                            f3      = fi[k2]    - a;
                            f2      = fi[k2]    + a;
                            g3      = gi[k2]    - b;
                            g2      = gi[k2]    + b;
                            b       = s1*f2     - c1*g3;
                            a       = c1*f2     + s1*g3;
                            fi[k2]  = f0        - a;
                            fi[0 ]  = f0        + a;
                            gi[k3]  = g1        - b;
                            gi[k1]  = g1        + b;
                            b       = c1*g2     - s1*f3;
                            a       = s1*g2     + c1*f3;
                            gi[k2]  = g0        - a;
                            gi[0 ]  = g0        + a;
                            fi[k3]  = f1        - b;
                            fi[k1]  = f1        + b;
                            gi     += k4;
                            fi     += k4;
                        } while (fi<fn);
                    }
                    TRIG_RESET(k,c1,s1);
                } while (k4<n);
            }
            
            void ifft64(int n, double *real, double *imag)
            {
                double a,b,c,d;
                double q,r,s,t;
                int i,j,k;
                fht64(real,n);
                fht64(imag,n);
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i]; b = real[j];  q=a+b; r=a-b;
                    c = imag[i]; d = imag[j];  s=c+d; t=c-d;
                    imag[i] = (s+r)*0.5;  imag[j] = (s-r)*0.5;
                    real[i] = (q-t)*0.5;  real[j] = (q+t)*0.5;
                }
            }
            
            void realfft64(int n, double *real)
            {
                double a,b;
                int i,j,k;
                
                fht64(real,n);
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i];
                    b = real[j];
                    real[j] = (a-b)*0.5;
                    real[i] = (a+b)*0.5;
                }
            }
            
            void fft64(int n, double *real, double *imag)
            {
                double a,b,c,d;
                double q,r,s,t;
                int i,j,k;
                
                #ifdef DBG
                debug_printf("fft real %.2f %.2f",
                             real[0],real[1]);
                
                debug_printf("fft imag %.2f %.2f",
                             imag[0],imag[1]);
                #endif
                
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i]; b = real[j];  q=a+b; r=a-b;
                    c = imag[i]; d = imag[j];  s=c+d; t=c-d;
                    real[i] = (q+t)*.5; real[j] = (q-t)*.5;
                    imag[i] = (s-r)*.5; imag[j] = (s+r)*.5;
                }
                fht64(real,n);
                fht64(imag,n);
                #ifdef DBG
                debug_printf("fft after real %.2f %.2f",
                             real[0],real[1]);
                
                debug_printf("fft after imag %.2f %.2f",
                             imag[0],imag[1]);
                #endif
            }
            
            void realifft64(int n, double *real)
            {
                double a,b;
                int i,j,k;
                for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
                    a = real[i];
                    b = real[j];
                    real[j] = (a-b);
                    real[i] = (a+b);
                }
                fht64(real,n);
            }
        };
    };
    
    ///****************************************************************************
    
    namespace fft {
        
        /// Simplest FFT class with fft function, ifft function, and modular bin number.
        /// Default number of bins = 2048
        
        /// @tparam T     render input as this datatype. algorithm was designed to assume the use of floating point.
        
        template<typename T = double>
        class simpleFFT : public fftSupport::fftTransforms {
        public:
            
            /// simpleFFT constructor and destructor
            /// @param fftsize     Size of FFT (number of bins)
            
            simpleFFT(int fftsize, bool inverse)
            {
                f_fftsize = fftsize;
                isInverse = inverse;
            }
            
            ~simpleFFT()
            {
                
            }
            
        
            /// DSP Setup message
            /// Pass this through the message<> dspsetup function of your object
        
            void dspsetup() {
                
                f_realin = 0L; // init pointers just to be safe
                f_imagin = 0L;
                f_realout = 0L;
                f_imagout = 0L;
                f_realin = (double *)max::sysmem_newptr(f_fftsize * 4 * sizeof(double));
                f_imagin = f_realin + f_fftsize;
                f_realout = f_realin + (f_fftsize*2);
                f_imagout = f_realin + (f_fftsize*3);
                
                f_interval = f_fftsize;

                f_countdown = 0;
                f_realinptr = f_realin;
                f_imaginptr = f_imagin;
                    
                // zero buffers in case they have junk in them (i.e. now makes no glitches when dsp is turned on)
                memset(f_realin, 0, sizeof(*f_realin) * f_fftsize*4);
                
            };

            
            
            /// Calculate FFT or IFFT every time the audio thread is called
            /// Extension of the object's operator()() method
            /// @param T sample  Pass samples through your object's operator()()
            ///                  function one by one
            
            void operator()(T sample) {

                //TODO: check if bin number argument changed. If it's changed,
                //replace circular buffer with another circular buffer equal to
                //the new size of the fft
                
                if (sampleCount < f_fftsize) {
                    f_buffer[sampleCount] = sample;
                    sampleCount++;
                } else {
                    //TODO: IFFT
                    if (isInverse) {
                        
                        memcpy(f_imaginptr, &f_buffer, sizeof(f_buffer));
                        memcpy(f_imagout, f_imagin, sizeof(*f_imagin) * f_fftsize*2);
      
                        ///This is where the magic happens
                        ifft((int) f_fftsize, (float*) f_realout, (float*) f_imagout);
                        
                    } else {
                                
                        memcpy(f_realinptr, &f_buffer, sizeof(f_buffer));
                        memcpy(f_realout, f_realin, sizeof(*f_realin) * f_fftsize*2);
                        
                        ///This is where the magic happens
                        realfft64((int) f_fftsize, f_realout);
                    }

                    // reset pointers
                    f_realinptr = f_realin;
                    f_imaginptr = f_imagin;
                    
                    // reset sampleCount
                    sampleCount = 0;
                    
                }
            };
            
            
            
            /// Calculate FFT or IFFT every time the audio thread is called
            /// Extension of the object's operator()() method
            /// @return       Post-transformed audio buffer in the form
            ///               of a buffer~ reference
            
            T magnitude(int index){
                return (sqrt((f_imagout[index] * f_imagout[index]) + (f_realout[index] * f_realout[index])));
            }
            
            
            
            /// Set the FFT size (i.e. number of bins)
            /// @param size       Set FFT size. FFT size must be between
            ///                   128 and 2048 and be a power of 2
            
            void setSize(int size){
                if(size > 128 && size < 2048){
                    if((size & (size - 1)) == 0){
                        f_fftsize = size;
                    } else {
                        error("Error: FFT size must be a power of 2");
                    }
                } else {
                    error("Error: FFT size must be between 128 and 2048");
                }
            }
            
            
            
            /// Get the FFT size (i.e. number of bins)
            /// @return       FFT size
            
            int getSize(){
                return f_fftsize;
            }
            
        
        private:
            bool isInverse = false;
            int	 f_fftsize = 2048;		// size
            int  sampleCount = 0;       // incrementor for f_buffer[] when audio is on
            
            long f_interval;		// sample count before doing another FFT (i.e. hop size, but always >= to fftsize)
            long f_countdown;	// vector count for phase offset fft
            double f_1overpts;		// (1.0/fftsize) for inverse FFT scaling
            
            double f_buffer[2048];  // audio buffer
            
            int m_minBins = 128;
            int m_maxBins = 2048;
            
            double *f_realin;		// where transform is done
            double *f_imagin;		// where transform is done
            double *f_realout;		// where transform is done
            double *f_imagout;		// where transform is done
            
            double *f_realinptr;	// input fill ptr into realin
            double *f_imaginptr;	// input fill ptr into imagin

        };

        
    
        /// Extended FFT class with more options including window types, transforms, and use of multiple
        /// FFTs in a single instance
        /// @tparam T       render output as this datatype. algorithm was designed to assume the use of floating point.
        
        template<typename T = buffer_lock<>>
        class extendedFFT : public fftSupport::fftTransforms {
        public:
        
        
        
        private:
        
        };
        
    };
    
    
    
}}}    // namespace c74::min::lib
