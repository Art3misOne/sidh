
package sidh;

/**************************************************************************************************
 *
 * Implements classes for Montgomery curves and points on the curve over GF(p) and GF(p^2).
 * There are three classes for points on the curve: over Fp in both projective and affine
 * coordinates and over Fp^2 in projective coordinates.
 *  
 **************************************************************************************************/

import java.math.BigInteger;

class FpPointProjective {
  /* Points with coordinates in Fp in projective coordinates */

  private final Felm x;
  private final Felm z;

  public static final FpPointProjective INFINITY = new FpPointProjective (Felm.ONE, Felm.ZERO);


  public FpPointProjective (Felm xc, Felm zc) {
    x = xc;
    z = zc;
  }


  public FpPointProjective (BigInteger xc, BigInteger zc) {
    x = new Felm (xc);
    z = new Felm (zc);
  }


  public FpPointProjective (FpPointProjective p) {
    x = p.x;
    z = p.z;
  }


  public Felm getX() {
    return x;
  }


  public Felm getZ() {
    return z;
  }


  public FpPointProjective swapPointsBasefield (FpPointProjective q, BigInteger option) {
    // If option = 0 then this <- this and q <- q, else this <- q and q <- this
    Felm qx, qz;

    qx = x.fpSwap(q.x, option);
    qz = z.fpSwap(q.z, option);

    return new FpPointProjective (qx, qz);
  }


  public String toString() {
    return "(" + x + ", " + z + ")";
  }
}


class FpPointAffine {
  /* Points with coordinates in Fp in affine coordinates */

  private final Felm x;
  private final Felm y;


  public FpPointAffine (Felm xc, Felm yc) {
    x = xc;
    y = yc;
  }

  public FpPointAffine (BigInteger xc, BigInteger yc) {
    x = new Felm (xc);
    y = new Felm (yc);
  }


  public FpPointAffine (FpPointAffine p) {
    x = p.x;
    y = p.y;
  }


  public Felm getX() {
    return x;
  }


  public Felm getY() {
    return y;
  }


  public FpPointAffine swapPointsBasefield (FpPointAffine q, BigInteger option) {
    // If option = 0 then this <- this and q <- q, else this <- q and q <- this
    Felm qx, qy;

    qx = x.fpSwap(q.x, option);
    qy = y.fpSwap(q.y, option);

    return new FpPointAffine (qx, qy);
  }


  public String toString() {
    return "(" + x + ", " + y + ")";
  }
}


class F2Point {
  /* Points with coordinates in Fp^2 in projective coordinates */

  private final F2elm x;
  private final F2elm z;

  public static final F2Point INFINITY = new F2Point (F2elm.ONE, F2elm.ZERO);


  public F2Point (F2elm xc, F2elm zc) {
    x = xc;
    z = zc;
  }


  public F2Point (Felm x0, Felm x1, Felm z0, Felm z1) {
    x = new F2elm (x0, x1);
    z = new F2elm (z0, z1);
  }

    
  public F2Point (F2Point p) {
    x = p.x;
    z = p.z;
  }


  public F2elm getX() {
    return x;
  }


  public F2elm getZ() {
    return z;
  }


  public F2Point swapPoints (F2Point q, BigInteger option) {
    // If option = 0 then this <- this and q <- q, else this <- q and q <- this
    F2elm qx, qz;

    qx = x.f2Swap(q.x, option);
    qz = z.f2Swap(q.z, option);

    return new F2Point (qx, qz);
  }


  public String toString() {
    return x.f2Mult(z.f2Inverse()).toString();
    //return "(" + x + ", " + z + ")";
  }
}


class MontCurve {
  /* 
   * Montgomery curves are of the form By^2 = Cx^3 + Ax^2 + Cx (using projective coefficients). 
   * Our computations do no require the parameter B, but it has been left in for completeness. 
   */

  protected F2elm a;
  protected F2elm b;
  protected F2elm c;

  // (a24/c24) = [(a/c) + 2]/4. Used often so it's advantageous to have precomputed.
  protected F2elm a24;              
  protected F2elm c24;


  public MontCurve() {
    // By default set a = 0, b = c = 1. Curve: y^2 = x^3 + x. a24 = 1 and c24 = 2.

    a = F2elm.ZERO;
    b = F2elm.ONE;
    c = F2elm.ONE;

    a24 = F2elm.ONE;
    c24 = c.f2Add (c);
  }


  public MontCurve (F2elm aA, F2elm cC) {
    a = aA;
    b = F2elm.ONE;
    c = cC;

    updateA24C24();
  }


  public MontCurve (F2elm aA, F2elm bB, F2elm cC) {
    a = aA;
    b = bB;
    c = cC;

    updateA24C24();
  }


  public MontCurve (MontCurve curveIn) {
    a = curveIn.a;
    b = curveIn.b;
    c = curveIn.c;
    a24 = curveIn.a24;
    c24 = curveIn.c24;
  }


  public void updateA24C24() {
    if (a.f2IsEven()) {
      a24 = c.f2Add (a.f2Div2());           // a24 = c + a/2
      c24 = c.f2Add (c);                     // c24 = 2c
    }
    else {
      c24 = c.f2Add (c);                     // c24 = 2c
      a24 = c24.f2Add (a);                   // a24 = 2c + a
      c24 = c24.f2Add (c24);                 // c24 = 4c
    }
  }


  public void normalizeA24C24() {
    a24 = a24.f2Mult(c24.f2Inverse());
    c24 = F2elm.ONE;
  }


  public F2elm getA() {
    return a;
  }


  public F2elm getC() {
    return c;
  }


  public static F2Point distortAndDiff (Felm px) {
    // Compute projective coordinates of difference point Q-P where Q = tau(P).
    // Inputs: px, the x-coordinate of an affine point
    // Outputs: projective coordinates (xd, zd) of Q-P 

    Felm xd1, zd0;

    xd1 = px.fpSqr();                    // xd = i*(px^2)
    xd1 = xd1.fpAdd(Felm.ONE);           // xd = i*(1+px^2)
    zd0 = px.fpAdd(px);                  // zd = px + px
    
    return new F2Point (Felm.ZERO, xd1, zd0, Felm.ZERO);
  }


  public FpPointProjective[] xDblAddBasefield (FpPointProjective p, FpPointProjective q, Felm xpq) {
    // Double and add in the base field. 
    // Inputs: Projective Fp points p (this) and q, affine difference xpq
    // Outputs: an array of two points, 2*p and p+q

    Felm t0, t1, t2, qx, qz, px, pz, fpC24, fpA24;
    FpPointProjective res[] = new FpPointProjective[2];

    fpC24 = c24.f2Get0();
    fpA24 = a24.f2Get0();

    qx = q.getX();
    qz = q.getZ();

    px = p.getX();
    pz = p.getZ();

    t0 = px.fpAdd(pz);                      // t0 = px+pz
    t1 = px.fpSub(pz);                      // t1 = px-pz
    px = t0.fpSqr();                        // px = (px+pz)^2
    t2 = qx.fpSub(qz);                      // t2 = qx-qz
    qx = qx.fpAdd(qz);                      // qx = qx+qz
    t0 = t0.fpMult(t2);                     // t0 = (px+pz)*(qx-qz)
    pz = t1.fpSqr();                        // pz = (px-pz)^2
    t1 = t1.fpMult(qx);                     // t1 = (px-pz)*(qx+qz)
    t2 = px.fpSub(pz);                      // t2 = (px+pz)^2-(px-pz)^2
    t2 = t2.fpMult(fpA24);                  // t2 = a24 * [(px+pz)^2-(px-pz)^2]

    pz = pz.fpMult(fpC24);                  // pz = c24*(px-pz)^2
    px = px.fpMult(pz);                     // px = c24*(px+pz)^2*(px-pz)^2
    pz = t2.fpAdd(pz);                      // pz = a24*[(px+pz)^2-(px-pz)^2]+c24*(px-pz)^2
    
    qz = t0.fpSub(t1);                      // qz = (px+pz)*(qx-qz)-(px-pz)*(qx+qz)
    qx = t0.fpAdd(t1);                      // qx = (px+pz)*(qx-qz)+(px-pz)*(qx+qz)
    pz = pz.fpMult(t2);                     // pz = [a24*[(px+pz)^2-(px-pz)^2]+c24*(px-pz)^2]
                                                //      * [(px+pz)^2-(px-pz)^2]
    qz = qz.fpSqr();                        // qz = [(px+pz)*(qx-qz)-(px-pz)*(qx+qz)]^2
    qx = qx.fpSqr();                        // qx = [(px+pz)*(qx-qz)+(px-pz)*(qx+qz)]^2
    qz = qz.fpMult(xpq);                    // qz = xpq*[(px+pz)*(qx-qz)-(px-pz)*(qx+qz)]^2    

    res[0] = new FpPointProjective (px, pz);
    res[1] = new FpPointProjective (qx, qz);

    return res;
  }


  public static F2Point xAdd (F2Point p, F2Point q, F2elm xpq) {
    // Differential addition of p and q 
    // Inputs: points p and q, and the xcoord of their affine difference xpq
    // Outputs: return p + q

    F2elm t0, t1, px, pz, qx, qz;

    px = p.getX();
    pz = p.getZ();

    qx = q.getX();
    qz = q.getZ();

    t0 = px.f2Add(pz);                      // t0 = px + pz
    t1 = px.f2Sub(pz);                      // t1 = px - pz
    px = qx.f2Sub(qz);                      // px = xq - qz
    pz = qx.f2Add(qz);                      // pz = xq + qz
    t0 = t0.f2Mult(px);                     // t0 = (px+pz)*(xq-qz)
    t1 = t1.f2Mult(pz);                     // t1 = (px-pz)*(xq+qz)
    pz = t0.f2Sub(t1);                      // pz = (px+pz)*(xq-qz) - (px-pz)*(xq+qz)
    px = t0.f2Add(t1);                      // px = (px+pz)*(xq-qz) - (px-pz)*(xq+qz)
    pz = pz.f2Sqr();                        // pz = [(px+pz)*(xq-qz) - (px-pz)*(xq+qz)]^2
    px = px.f2Sqr();                        // px = [(px+pz)*(xq-qz) - (px-pz)*(xq+qz)]^2
    pz = pz.f2Mult(xpq);                    // pz = pxq*[(px+pz)*(xq-qz)-(px-pz)*(xq+qz)]^2

    return new F2Point (px, pz);
  }


  public F2Point xDbl (F2Point p) {
    // Double the point p
    // Input: Projective point p with coordinates in Fp^2
    // Output: (qx:qz) = 2*(px:pz)

    F2elm t0, t1, qx, qz, px, pz;

    px = p.getX();
    pz = p.getZ();

    t0 = px.f2Sub(pz);                   // t0 = px - pz
    t1 = px.f2Add(pz);                   // t1 = px + pz
    t0 = t0.f2Sqr();                     // t0 = (px-pz)^2
    t1 = t1.f2Sqr();                     // t1 = (px+pz)^2
    qz = c24.f2Mult(t0);                 // qz = c24*(px-pz)^2
    qx = t1.f2Mult(qz);                  // qx = c24*(px-pz)^2*(px+pz)^2
    t1 = t1.f2Sub(t0);                   // t1 = (px+pz)^2 - (px-pz)^2
    t0 = a24.f2Mult(t1);                 // t0 = a24 * [(px+pz)^2 - (px-pz)^2]
    qz = qz.f2Add(t0);                   // qz = a24*[(px+pz)^2-(px-pz)^2] + c24*(px-pz)^2
    qz = qz.f2Mult(t1);                  // qz = [a24*[(px+pz)^2-(px-pz)^2] + c24*(px-pz)^2]
                                         //      * [(px+pz)^2-(px-pz)^2]

    return new F2Point (qx, qz);
  }


  public F2Point xDble (F2Point p, int e) {
    // Computes [2^e](px:pz) via e repeated doublings
    // Inputs: F2Point p and exponent e 
    // Output: (qx:qz) = (2^e)*(px:pz)

    F2Point q;
    int i;

    q = p;

    for (i = 0; i < e; i++)
      q = xDbl(q);

    return q;
  }


  public F2Point xTpl (F2Point p) {
    // Given point p compute 3*p by using 3p = 2p + p and doing efficient double and add
    // Inputs: F2Point p
    // Outputs: (qx:qz) = 3*(px:pz)

    F2elm t0, t1, t2, t3, t4, t5, x, z;

    x = p.getX();
    z = p.getZ();

    t2 = x.f2Add(z);                     // t2 = x+z
    t0 = x.f2Sqr();                      // t0 = x^2
    t1 = z.f2Sqr();                      // t1 = z^2
    t2 = t2.f2Sqr();                     // t2 = (x+z)^2
    t3 = t0.f2Mult(c);                   // t3 = c*x^2
    t2 = t2.f2Sub(t0);                   // t2 = (x+z)^2-x^2
    t4 = t1.f2Mult(c);                   // t4 = c*z^2
    t2 = t2.f2Sub(t1);                   // t2 = 2xz = (x+z)^2-x^2-z^2
    t5 = t3.f2Add(t4);                   // t5 = c*x^2+c*z^2
    t2 = t2.f2Mult(a);                   // t2 = 2axz
    t3 = t3.f2Add(t3);                   // t3 = 2c*x^2 
    t4 = t4.f2Add(t4);                   // t4 = 2c*z^2
    t3 = t3.f2Add(t2);                   // t3 = 2c*x^2+2axz
    t4 = t4.f2Add(t2);                   // t4 = 2c*z^2+2axz
    t3 = t3.f2Add(t5);                   // t3 = 2c*x^2+2axz + c*x^2+c*z^2
    t4 = t4.f2Add(t5);                   // t4 = 2c*z^2+2axz + c*x^2+c*z^2
    t2 = t0.f2Sub(t1);                   // t2 = x^2-z^2
    t0 = t0.f2Add(t0);                   // t0 = 2x^2
    t1 = t1.f2Add(t1);                   // t1 = 2z^2
    t2 = t2.f2Mult(t5);                  // t2 = (x^2-z^2)(c*x^2+c*z^2)
    t1 = t1.f2Mult(t3);                  // t1 = 2z^2*[2c*x^2+2axz+c*x^2+c*z^2]
    t0 = t0.f2Mult(t4);                  // t0 = 2x^2*[2c*z^2+2axz+c*x^2+c*z^2]
    t1 = t1.f2Sub(t2);                   // t1 = 2z^2*[2c*x^2+2axz+c*x^2+c*z^2] 
                                         //      - (x^2-z^2)(c*x^2+c*z^2)
    t0 = t0.f2Add(t2);                   // t0 = 2x^2*[2c*z^2+2axz+c*x^2+c*z^2] 
                                         //      + (x^2-z^2)(c*x^2+c*z^2)
    t1 = t1.f2Sqr();                     // t1 = [2z^2*[2c*x^2+2axz+c*x^2+c*z^2] 
                                         //      - (x^2-z^2)(c*x^2+c*z^2)]^2
    t0 = t0.f2Sqr();                     // t0 = [2x^2*[2c*z^2+2axz+c*x^2+c*z^2] 
                                         //      + (x^2-z^2)(c*x^2+c*z^2)]^2
    t3 = x.f2Mult(t1);                   // t3 = x*[2z^2*[2c*x^2+2axz+c*x^2+c*z^2] 
                                         //      - (x^2-z^2)(c*x^2+c*z^2)]^2
    t4 = z.f2Mult(t0);                   // t4 = z*[2x^2*[2c*z^2+2axz+c*x^2+c*z^2] 
                                         //      + (x^2-z^2)(c*x^2+c*z^2)]^2

    return new F2Point (t3, t4);
  }


  public F2Point xTple (F2Point p, int e) {
    // Computes [3^e](px:pz) via e repeated triplings
    // Input: curve point P, exponent e
    // Output: q = [3^e]*P 
    
    int i;

    F2Point q = p;

    for (i = 0; i < e; i++)
      q = xTpl (q);

    return q;
  }


  public F2Point[] xDblAdd (F2Point p, F2Point q, F2elm xpq) {
    // Simultaneous double and differential addition.
    // Inputs: projective points p and q and the x-coordinate of their affine difference xpq
    // Outputs: 2*p and p+q

    F2elm t0, t1, t2, t3, qx, qz, px, pz;
    F2Point pq[];

    px = p.getX();
    pz = p.getZ();
    qx = q.getX();
    qz = q.getZ();

    t0 = px.f2Add(pz);                   // t0 = px+pz
    t1 = px.f2Sub(pz);                   // t1 = px-pz
    px = t0.f2Sqr();                     // px = (px+pz)^2
    t2 = qx.f2Sub(qz);                   // t2 = qx-qz
    qx = qx.f2Add(qz);                   // qx = qx+qz
    t0 = t0.f2Mult(t2);                  // t0 = (px+pz)*(qx-qz)
    pz = t1.f2Sqr();                     // pz = (px-pz)^2
    t1 = t1.f2Mult(qx);                  // t1 = (px-pz)*(qx+qz)
    t3 = px.f2Sub(pz);                   // t2 = (px+pz)^2 - (px-pz)^2

    t2 = t3.f2Mult(a24);                 // t2 = a24 * [(px+pz)^2-(px-pz)^2]
    pz = pz.f2Mult(c24);                 // pz = c24 * (px-pz)^2
    px = px.f2Mult(pz);                  // px = c24*(px-pz)^2 * (px+pz)^2
    pz = t2.f2Add(pz);                   // pz = a24*[(px+pz)^2-(px-pz)^2] + c24*(px-pz)^2

    qz = t0.f2Sub(t1);                   // qz = (px+pz)*(qx-qz) - (px-pz)*(qx+qz)
    qx = t0.f2Add(t1);                   // qx = (px+pz)*(qx-qz) + (px-pz)*(qx+qz)
    pz = pz.f2Mult(t3);                  // pz = [a24*[(px+pz)^2-(px-pz)^2]+c24*(px-pz)^2] 
                                             //      * [(px+pz)^2-(px-pz)^2]
    qz = qz.f2Sqr();                     // qz = [(px+pz)*(qx-qz)-(px-pz)*(qx+qz)]^2
    qx = qx.f2Sqr();                     // qx = [(px+pz)*(qx-qz)+(px-pz)*(qx+qz)]^2
    qz = qz.f2Mult(xpq);                 // qz = xpq * [(px+pz)*(qx-qz)-(px-pz)*(qx+qz)]^2

    pq = new F2Point[2];
    pq[0] = new F2Point (px, pz);
    pq[1] = new F2Point (qx, qz);

    return pq;
  }


  public FpPointProjective[] ladder (Felm x, BigInteger m, int obits) {
    // Montgomery ladder to compute a multiple of a point on a curve
    // Inputs: Affine x-coordinate of a point, scalar m
    // Outputs: m * (x:1)
    // Note: this method simultaneously computes m*(x:1) and (m+1)*(x:1) and there are some 
    //       computations which require both, so both are returned

    int nbits, i;
    BigInteger t, bit = BigInteger.ZERO;
    FpPointProjective pq[];
    
    pq = new FpPointProjective[2];

    // Initialize p and q with (1:0) and (x:1)
    pq[0] = new FpPointProjective (FpPointProjective.INFINITY);
    pq[1] = new FpPointProjective (x, Felm.ONE);

    nbits = m.bitLength();

    // Double/Add the identity up to the maximum order bits so that this method takes constant
    // time regardless of how many bits are in the scalar m (prevents timing attacks). 
    for (i = obits; i > nbits; i--) {
      pq[1] = pq[0].swapPointsBasefield (pq[1], bit);
      pq = xDblAddBasefield (pq[0], pq[1], x);
      pq[1] = pq[0].swapPointsBasefield (pq[1], bit);
    }

    // Actual computation: double/add for the bits in the scalar m
    for (i = nbits-1; i >= 0; i--) {
      bit = m.testBit(i) ? BigInteger.ONE : BigInteger.ZERO;

      // if bit = 0 then p = 2p and q = p+q
      // else if bit = 1 then q = 2q and p = p+q 
      pq[1] = pq[0].swapPointsBasefield (pq[1], bit);
      pq = xDblAddBasefield (pq[0], pq[1], x);
      pq[1] = pq[0].swapPointsBasefield (pq[1], bit);
    }

    return pq;
  }


  public F2Point ladder3pt (F2elm px, F2elm qx, F2elm xpq, BigInteger m, int obits) {
    // Computes P + m[Q] via x-only arithmetic.
    // Inputs: 3 affine x coordinates px, qx, xpq and scalar m
    // Outputs: (wx:wz) = P + m[Q] 

    F2Point w, uv[];
    F2elm temp1, temp2;
    int i, nbits;
    BigInteger bit = BigInteger.ZERO;

    uv = new F2Point[2];
    uv[0] = F2Point.INFINITY;                             // u = (1:0)
    uv[1] = new F2Point (qx, F2elm.ONE);                  // v = (qx:1)
    w = new F2Point (px, F2elm.ONE);                      // w = (px:1)

    nbits = m.bitLength();

    // Double/Add the identity up to the maximum order bits so that this method takes constant
    // time regardless of how many bits are in the scalar m (prevents timing attacks). 
    for (i = obits; i > nbits; i--) {
      uv[0] = w.swapPoints(uv[0], bit);
      uv[1] = uv[0].swapPoints(uv[1], bit);
      temp1 = px.f2Select(qx, bit);
      temp2 = qx.f2Select(xpq, bit);
      w = xAdd(w, uv[0], temp1);
      uv = xDblAdd(uv[0], uv[1], temp2); 
      uv[1] = uv[0].swapPoints(uv[1], bit);
      uv[0] = w.swapPoints(uv[0], bit);
    }

    for (i = nbits-1; i >= 0; i--) {
      bit = m.testBit(i) ? BigInteger.ONE : BigInteger.ZERO;

      // if bit = 0 then w = w+u, u = 2*u, and v = u+v
      // else if bit = 1 then u = u+v, v = 2*v, and w = v+w
      uv[0] = w.swapPoints(uv[0], bit);
      uv[1] = uv[0].swapPoints(uv[1], bit);
      temp1 = px.f2Select(qx, bit);
      temp2 = qx.f2Select(xpq, bit);
      w = xAdd(w, uv[0], temp1);
      uv = xDblAdd(uv[0], uv[1], temp2); 
      uv[1] = uv[0].swapPoints(uv[1], bit);
      uv[0] = w.swapPoints(uv[0], bit);
    }

    return w;
  }


  public F2elm jInv() {
    // Computes the j-invariant of a Montgomery curve with projective constant a/c. 
    // Curve: B*y^2=x^3+(a/c)*x^2+x
    // Outputs: j = 256*(a^2-3*c^2)^3/(c^4*(a^2-4*c^2))

    F2elm t0, t1, jinv;

    jinv = a.f2Sqr();                       // jinv = a^2
    t1 = c.f2Sqr();                         // t1 = c^2
    t0 = t1.f2Add(t1);                      // t0 = t1 + t1
    t0 = jinv.f2Sub(t0);                    // t0 = jinv - t0
    t0 = t0.f2Sub(t1);                      // t0 = t0 - t1
    jinv = t0.f2Sub(t1);                    // jinv = t0 - t1
    t1 = t1.f2Sqr();                        // t1 = t1^2
    jinv = jinv.f2Mult(t1);                 // jinv = jinv*t1
    t0 = t0.f2Add(t0);                      // t0 = t0 + t0
    t0 = t0.f2Add(t0);                      // t0 = t0 + t0
    t1 = t0.f2Sqr();                        // t1 = t0^2
    t0 = t0.f2Mult(t1);                     // t0 = t0 * t1
    t0 = t0.f2Add(t0);                      // t0 = t0 + t0
    t0 = t0.f2Add(t0);                      // t0 = t0 + t0
    jinv = jinv.f2Inverse();                // jinv = 1/jinv
    jinv = t0.f2Mult(jinv);                 // jinv = t0*jinv

    return jinv;
  }


  public String toString() {
    return b + " y^2 = " + c + " x^3 + " + a + " x^2 + " + c + " x   where (a24/c24) = (" + a24 + "/" + ")";
  }
}
