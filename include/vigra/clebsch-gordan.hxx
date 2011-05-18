/************************************************************************/
/*                                                                      */
/*        Copyright 2009-2010 by Ullrich Koethe and Janis Fehr          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_CLEBSCH_GORDAN_HXX
#define VIGRA_CLEBSCH_GORDAN_HXX

#include "config.hxx"
#include "numerictraits.hxx"
#include "error.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"

namespace vigra {

namespace {

void ThreeJSymbolM(double l1, double l2, double l3, double m1,
                   double &m2min, double &m2max, double *thrcof, int ndim,
                   int &errflag)
{
    ContractViolation err;  
    const double zero = 0.0, eps = 0.01, one = 1.0, two = 2.0;

    int nfin, nlim, i, n, index, lstep, nfinp1, nfinp2, nfinp3, nstep2;
    double oldfac, dv, newfac, sumbac = 0.0, thresh, a1s, sumfor, sumuni,
    sum1, sum2, x, y, m2, m3, x1, x2, x3, y1, y2, y3, cnorm,
    ratio, a1, c1, c2, c1old = 0.0, sign1, sign2;

    // Parameter adjustments
    --thrcof;

    errflag = 0;

    // "hugedouble" is the square root of one twentieth of the largest floating
    // point number, approximately.
    double hugedouble   = std::sqrt(NumericTraits<double>::max() / 20.0),
    srhuge = std::sqrt(hugedouble),
    tiny   = one / hugedouble,
    srtiny = one / srhuge;

    // lmatch = zero

    //  Check error conditions 1, 2, and 3.
    if (l1 - abs(m1) + eps < zero
            || std::fmod(l1 + abs(m1) + eps, one) >= eps + eps) 
    {
        errflag = 1;
        err << "ThreeJSymbolM: l1-abs(m1) less than zero or l1+abs(m1) not integer.\n";
        throw err;
    } 
    else if (l1+l2-l3 < -eps || l1-l2+l3 < -eps || -(l1) + l2+l3 < -eps) 
    {
        errflag = 2;
        err << " ThreeJSymbolM: l1, l2, l3 do not satisfy triangular condition:"
            << l1 << " " << l2 << " " << l3 << "\n";
        throw err;
    } 
    else if (std::fmod(l1 + l2 + l3 + eps, one) >= eps + eps) 
    {
        errflag = 3;
        err << " ThreeJSymbolM: l1+l2+l3 not integer.\n";
        throw err;
    }

    // limits for m2
    m2min = std::max(-l2,-l3-m1);
    m2max = std::min(l2,l3-m1);

    // Check error condition 4.
    if (std::fmod(m2max - m2min + eps, one) >= eps + eps) {
        errflag = 4;
        err << " ThreeJSymbolM: m2max-m2min not integer.\n";
        throw err;
    }
    if (m2min < m2max - eps) 
        goto L20;
    if (m2min < m2max + eps) 
        goto L10;

    //  Check error condition 5.
    errflag = 5;
    err << " ThreeJSymbolM: m2min greater than m2max.\n";
    throw err;

    // This is reached in case that m2 and m3 can take only one value.
L10:
    // mscale = 0
    thrcof[1] = (odd(int(abs(l2-l3-m1)+eps)) 
                       ? -one 
                       :  one) / std::sqrt(l1+l2+l3+one);
    return;

    // This is reached in case that M1 and M2 take more than one value.
L20:
    // mscale = 0
    nfin = int(m2max - m2min + one + eps);
    if (ndim - nfin >= 0) 
        goto L23;

    // Check error condition 6.

    errflag = 6;
    err << " ThreeJSymbolM: Dimension of result array for 3j coefficients too small.\n";
    throw err;

    //  Start of forward recursion from m2 = m2min

L23:
    m2 = m2min;
    thrcof[1] = srtiny;
    newfac = 0.0;
    c1 = 0.0;
    sum1 = tiny;

    lstep = 1;
L30:
    ++lstep;
    m2 += one;
    m3 = -m1 - m2;

    oldfac = newfac;
    a1 = (l2 - m2 + one) * (l2 + m2) * (l3 + m3 + one) * (l3 - m3);
    newfac = std::sqrt(a1);

    dv = (l1+l2+l3+one) * (l2+l3-l1) - (l2-m2+one) * (l3+m3+one)
         - (l2+m2-one) * (l3-m3-one);

    if (lstep - 2 > 0) 
        c1old = abs(c1);

    // L32:
    c1 = -dv / newfac;

    if (lstep > 2) 
        goto L60;

    //  If m2 = m2min + 1, the third term in the recursion equation vanishes,  
    //  hence

    x = srtiny * c1;
    thrcof[2] = x;
    sum1 += tiny * c1 * c1;
    if (lstep == nfin) 
        goto L220;
    goto L30;

L60:
    c2 = -oldfac / newfac;

    // Recursion to the next 3j coefficient
    x = c1 * thrcof[lstep-1] + c2 * thrcof[lstep-2];
    thrcof[lstep] = x;
    sumfor = sum1;
    sum1 += x * x;
    if (lstep == nfin) 
        goto L100;

    // See if last unnormalized 3j coefficient exceeds srhuge
    if (abs(x) < srhuge) 
        goto L80;

    // This is reached if last 3j coefficient larger than srhuge,
    // so that the recursion series thrcof(1), ... , thrcof(lstep)
    // has to be rescaled to prevent overflow

    // mscale = mscale + 1
    for (i = 1; i <= lstep; ++i) 
    {
        if (abs(thrcof[i]) < srtiny) 
            thrcof[i] = zero;
        thrcof[i] /= srhuge;
    }
    sum1 /= hugedouble;
    sumfor /= hugedouble;
    x /= srhuge;

    // As long as abs(c1) is decreasing, the recursion proceeds towards
    // increasing 3j values and, hence, is numerically stable.  Once
    // an increase of abs(c1) is detected, the recursion direction is
    // reversed.

L80:
    if (c1old - abs(c1) > 0.0) 
        goto L30;

    //  Keep three 3j coefficients aroundi mmatch for comparison later
    //  with backward recursion values.

L100:
    // mmatch = m2 - 1
    nstep2 = nfin - lstep + 3;
    x1 = x;
    x2 = thrcof[lstep-1];
    x3 = thrcof[lstep-2];

    //  Starting backward recursion from m2max taking nstep2 steps, so
    //  that forwards and backwards recursion overlap at the three points
    //  m2 = mmatch+1, mmatch, mmatch-1.

    nfinp1 = nfin + 1;
    nfinp2 = nfin + 2;
    nfinp3 = nfin + 3;
    thrcof[nfin] = srtiny;
    sum2 = tiny;

    m2 = m2max + two;
    lstep = 1;
L110:
    ++lstep;
    m2 -= one;
    m3 = -m1 - m2;
    oldfac = newfac;
    a1s = (l2-m2+two) * (l2+m2-one) * (l3+m3+two) * (l3-m3-one);
    newfac = std::sqrt(a1s);
    dv = (l1+l2+l3+one) * (l2+l3-l1) - (l2-m2+one) * (l3+m3+one)
          - (l2+m2-one) * (l3-m3-one);
    c1 = -dv / newfac;
    if (lstep > 2) 
        goto L120;

    // if m2 = m2max + 1 the third term in the recursion equation vanishes

    y = srtiny * c1;
    thrcof[nfin - 1] = y;
    if (lstep == nstep2) 
        goto L200;
    sumbac = sum2;
    sum2 += y * y;
    goto L110;

L120:
    c2 = -oldfac / newfac;

    // Recursion to the next 3j coefficient

    y = c1 * thrcof[nfinp2 - lstep] + c2 * thrcof[nfinp3 - lstep];

    if (lstep == nstep2) 
        goto L200;

    thrcof[nfinp1 - lstep] = y;
    sumbac = sum2;
    sum2 += y * y;

    // See if last 3j coefficient exceeds SRHUGE

    if (abs(y) < srhuge) 
        goto L110;

    // This is reached if last 3j coefficient larger than srhuge,
    // so that the recursion series thrcof(nfin), ... , thrcof(nfin-lstep+1)   
    // has to be rescaled to prevent overflow.

    // mscale = mscale + 1
    for (i = 1; i <= lstep; ++i) 
    {
        index = nfin - i + 1;
        if (abs(thrcof[index]) < srtiny) 
            thrcof[index] = zero;
        thrcof[index] /= srhuge;
    }
    sum2 /= hugedouble;
    sumbac /= hugedouble;

    goto L110;

    //  The forward recursion 3j coefficients x1, x2, x3 are to be matched
    //  with the corresponding backward recursion values y1, y2, y3.

L200:
    y3 = y;
    y2 = thrcof[nfinp2-lstep];
    y1 = thrcof[nfinp3-lstep];

    //  Determine now ratio such that yi = ratio * xi  (i=1,2,3) holds
    //  with minimal error.

    ratio = (x1*y1 + x2*y2 + x3*y3) / (x1*x1 + x2*x2 + x3*x3);
    nlim = nfin - nstep2 + 1;

    if (abs(ratio) < one) 
        goto L211;
    for (n = 1; n <= nlim; ++n)
        thrcof[n] = ratio * thrcof[n];
    sumuni = ratio * ratio * sumfor + sumbac;
    goto L230;

L211:
    ++nlim;
    ratio = one / ratio;
    for (n = nlim; n <= nfin; ++n)
        thrcof[n] = ratio * thrcof[n];
    sumuni = sumfor + ratio * ratio * sumbac;
    goto L230;
L220:
    sumuni = sum1;

    // Normalize 3j coefficients

L230:
    cnorm = one / std::sqrt((l1+l1+one) * sumuni);

    // Sign convention for last 3j coefficient determines overall phase

    sign1 = sign(thrcof[nfin]);
    sign2 = odd(int(abs(l2-l3-m1)+eps)) 
                   ? -one 
                   : one;
    if (sign1 * sign2 <= 0.0) 
        goto L235;
    else 
        goto L236;

L235:
    cnorm = -cnorm;

L236:
    if (abs(cnorm) < one) 
        goto L250;

    for (n = 1; n <= nfin; ++n)
        thrcof[n] = cnorm * thrcof[n];
    return;

L250:
    thresh = tiny / abs(cnorm);
    for (n = 1; n <= nfin; ++n) 
    {
        if (abs(thrcof[n]) < thresh) 
            thrcof[n] = zero;
        thrcof[n] = cnorm * thrcof[n];
    }
}

} // anonymous namespace

inline 
double clebschGordan (double l1, double m1, double l2, double m2, double l3, double m3)
{
    const double err = 0.01;
    double CG = 0.0, m2min, m2max, *cofp;                               
    // static array for calculation of 3-j symbols 
    const int ncof = 100;
    double cof[ncof];                                            
    // reset error flag   
    int errflag = 0;
    ContractViolation Err;

    // Check for physical restriction.
    // All other restrictions are checked by the 3-j symbol routine.
    if ( abs(m1 + m2 - m3) > err) 
    {
        errflag = 7;
        Err << " clebschGordan: m1 + m2 - m3 is not zero.\n";
        throw Err;
    }                                                                   
    // calculate minimum storage size needed for ThreeJSymbolM()
    // if the dimension becomes negative the 3-j routine will capture it
    int njm = roundi(std::min(l2,l3-m1) - std::max(-l2,-l3-m1) + 1);            

    // allocate dynamic memory if necessary
    ArrayVector<double> cofa;
    if(njm > ncof)
    {
        cofa.resize(njm);
        cofp = cofa.begin();
    }
    else
    {
        cofp = cof;
    }

    // calculate series of 3-j symbols
    ThreeJSymbolM (l1,l2,l3,m1, m2min,m2max, cofp,njm, errflag);        
    // calculated Clebsch-Gordan coefficient
    if (! errflag)
    {
        CG = cofp[roundi(m2-m2min)] * (odd(roundi(l1-l2+m3)) ? -1 : 1) * std::sqrt(2*l3+1);
    }
    else
    {
        Err << " clebschGordan: 3jM-sym error.\n";
        throw Err;
    }
    return CG;                                                          
}

} // namespace vigra 

#endif // VIGRA_CLEBSCH_GORDAN_HXX

