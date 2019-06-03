
package sidh;

/**************************************************************************************************
 *
 * Implements classes for Montgomery curves and points on the curve over GF(p) and GF(p^2).
 * There are three classes for points on the curve: over Fp in both projective and affine
 * coordinates and over Fp^2 in projective coordinates.
 *  
 **************************************************************************************************/

import java.math.BigInteger;


class F2Point {
  /* Points with coordinates in Fp^2 using (x:z) coordinates */

  private F2elm x;
  private F2elm z;

  public static final F2Point INFINITY = new F2Point (F2elm.ONE, F2elm.ZERO);


  public F2Point (F2elm xc, F2elm zc) {
    x = new F2elm (xc);
    z = new F2elm (zc);
  }


  public F2Point (Felm x0, Felm x1, Felm z0, Felm z1) {
    x = new F2elm (x0, x1);
    z = new F2elm (z0, z1);
  }

    
  public F2Point (F2Point p) {
    x = new F2elm (p.x);
    z = new F2elm (p.z);
  }


  public F2elm getX() {
    return x;
  }


  public F2elm getZ() {
    return z;
  }


  public void setX (F2elm xval) {
    x = xval;
  }


  public F2Point swapPoints (F2Point q, BigInteger option) {
    // If option = 0 then this <- this and q <- q, else this <- q and q <- this
    F2elm qx, qz;

    qx = x.f2Swap(q.x, option);
    qz = z.f2Swap(q.z, option);

    return new F2Point (qx, qz);
  }


  public void normalize () {
    z.f2InverseInPlace ();
    x.f2MultInPlace (z);
    z = new F2elm (F2elm.ONE);
  }
    

  public String toString() {
    return "(" + x + ": " + z + ")";
  }
}


class MontCurve {
  /* 
   * Montgomery curves are of the form By^2 = x^3 + (A/C)x^2 + x (using projective coefficients).
   * It is assumed that B = 1 so it is omitted.
   */

  protected F2elm a;
  protected F2elm c;

  // Useful precomputed values:

  protected F2elm aMinus2c;
  protected F2elm aPlus2c;

  protected F2elm a24;              
  protected F2elm c4;              
        

  public MontCurve() {
    // Default curve: y^2 = x^3 + x, ie set a = 0, b = c = 1.
    // Precomputed values are set when needed

    a = new F2elm (F2elm.ZERO);
    c = new F2elm (F2elm.ONE);
  }


  public MontCurve (F2elm aA, F2elm cC) {
    a = aA;
    c = cC;
  }


  public MontCurve (MontCurve curveIn) {
    a = new F2elm (curveIn.a);
    c = new F2elm (curveIn.c);
  }


  public void initializeConstants () {
    updateA24();
    updateC4();
    updatePlusMinus();  
  }

    
  public void updateA24 () {
    a24 = F2elm.leftShift (F2elm.ONE, 1);
    a24.f2AddInPlace (a);
    a24.f2Div2InPlace ();
    a24.f2Div2InPlace ();
  }


  public void updateC4 () {
    c4 = F2elm.leftShift (c, 2);
  }


  public void updatePlusMinus () {
    F2elm c2 = F2elm.leftShift (c, 1);
    aPlus2c = F2elm.add (a, c2);
    aMinus2c = F2elm.sub (a, c2);
  }
    

  public void updateAC (int order) {
    if (order == 3) {
      a = F2elm.add (aPlus2c, aMinus2c);
      a.f2LeftShiftInPlace (1);
      c = F2elm.sub (aPlus2c, aMinus2c);
    }

    else {
      c = F2elm.div2 (c4);
      a = F2elm.sub (aPlus2c, c);
      c.f2Div2InPlace ();
    }
  }


  public F2elm getAPlus2c () {
    return aPlus2c;
  }

    
  public void setAPlus2c (F2elm invalue) {
    aPlus2c = invalue;
  }


  public void setAMinus2c (F2elm invalue) {
    aMinus2c = invalue;
  }


  public void setC4 (F2elm invalue) {
    c4 = invalue;
  }


  public static F2elm recoverA (F2elm px, F2elm qx, F2elm dx) {
    F2elm t0, t1, ra;

    t1 = F2elm.add (px, qx);
    t0 = F2elm.mult (px, qx);
    ra = F2elm.mult (dx, t1);
    ra.f2AddInPlace (t0);
    t0.f2MultInPlace (dx);
    ra.f2SubInPlace (F2elm.ONE);
    t0.f2LeftShiftInPlace (2);
    t1.f2AddInPlace (dx);
    ra.f2SqrInPlace ();
    t0.f2InverseInPlace ();
    ra.f2MultInPlace (t0);
    ra.f2SubInPlace (t1);

    return ra;
  }
    

  public F2Point xDbl (F2Point p) {
    F2elm t0, t1, qx, qz, px, pz;

    px = p.getX();
    pz = p.getZ();

    t0 = F2elm.sub (px, pz);
    t1 = F2elm.add (px, pz);
    t0.f2SqrInPlace ();
    t1.f2SqrInPlace ();
    qz = F2elm.mult (c4, t0);
    qx = F2elm.mult (t1, qz);
    t1.f2SubInPlace (t0);
    t0 = F2elm.mult (aPlus2c, t1);
    qz.f2AddInPlace (t0);
    qz.f2MultInPlace (t1);

    return new F2Point (qx, qz);
  }


  public F2Point xDble (F2Point p, int e) {
    // Computes [2^e](px:pz) via e repeated doublings

    F2Point q = p;
    int i;

    for (i = 0; i < e; i++)
      q = xDbl (q);

    return q;
  }


  public F2Point xTpl (F2Point p) {
    // Given point p compute 3*p using the point tripling algorithm in "A Faster Software
    // Implementation of the Supersingular Isogeny Diffie-Hellman Key Exchange Protocol" by
    // Faz-Hernandez, Lopez, Ochoa-Jimenez, Rodriguez-Henriquez

    F2elm t0, t1, t2, t3, t4, t5, t6, px, pz, qx, qz;

    px = p.getX();
    pz = p.getZ();

    t0 = F2elm.sub (px, pz);
    t2 = F2elm.sqr (t0);
    t1 = F2elm.add (px, pz);
    t3 = F2elm.sqr (t1);
    t4 = F2elm.leftShift (px, 1);
    t0 = F2elm.leftShift (pz, 1);
    t1 = F2elm.sqr (t4);
    t1.f2SubInPlace (t3);
    t1.f2SubInPlace (t2);
    t5 = F2elm.mult (t3, aPlus2c);
    t3.f2MultInPlace (t5);
    t6 = F2elm.mult (t2, aMinus2c);
    t2.f2MultInPlace (t6);
    t3 = F2elm.sub (t2, t3);
    t2 = F2elm.sub (t5, t6);
    t1.f2MultInPlace (t2);
    t2 = F2elm.add (t1, t3);
    t2.f2SqrInPlace ();
    qx = F2elm.mult (t4, t2);
    t1 = F2elm.sub (t3, t1);
    t1.f2SqrInPlace ();
    qz = F2elm.mult (t0, t1);
    
    return new F2Point (qx, qz);
  }


  public F2Point xTple (F2Point p, int e) {
    // Computes [3^e](px:pz) via e repeated triplings
    
    int i;

    F2Point q = p;

    for (i = 0; i < e; i++)
      q = xTpl (q);

    return q;
  }


  public F2Point[] xDblAdd (F2Point p, F2Point q, F2elm xpq) {
    // Simultaneous double and differential addition.
    // Outputs: Array of points = [2*p, p+q]

    F2elm t0, t1, t2, qx, qz, px, pz;
    F2Point pq[];

    px = new F2elm (p.getX());
    pz = new F2elm (p.getZ());
    qx = new F2elm (q.getX());
    qz = new F2elm (q.getZ());

    t0 = F2elm.add (px, pz);
    t1 = F2elm.sub (px, pz);
    px = F2elm.sqr (t0);    
    t2 = F2elm.sub (qx, qz);
    qx.f2AddInPlace (qz);   
    t0.f2MultInPlace (t2);  
    pz = F2elm.sqr (t1);    
    t1.f2MultInPlace (qx);  
    t2 = F2elm.sub (px, pz);
    px.f2MultInPlace (pz);
    qx = F2elm.mult (t2, a24);
    qz = F2elm.sub (t0, t1); 
    pz.f2AddInPlace (qx);              
    qx = F2elm.add (t0, t1);            
    pz.f2MultInPlace (t2);                            
    qz.f2SqrInPlace ();                 
    qx.f2SqrInPlace ();                 
    qz.f2MultInPlace (xpq);             

    pq = new F2Point[2];
    pq[0] = new F2Point (px, pz);
    pq[1] = new F2Point (qx, qz);

    return pq;
  }


  public F2Point ladder3pt (F2elm xp, F2elm xq, F2elm xpq, BigInteger m, int obits) {
    // Computes P + m[Q] via x-only arithmetic.
    
    F2Point rs[], r;
    int i;
    BigInteger swap, bit, prevbit = BigInteger.ZERO;
    F2elm xval;
    
    rs = new F2Point[2];

    rs[0] = new F2Point (xq, F2elm.ONE);
    rs[1] = new F2Point (xpq, F2elm.ONE);
    r = new F2Point (xp, F2elm.ONE);

    for (i = 0; i < obits; i++) {
      bit = m.testBit(i) ? BigInteger.ONE : BigInteger.ZERO;
      swap = bit.xor(prevbit);
      prevbit = bit;

      r = rs[1].swapPoints (r, swap);
      rs = xDblAdd(rs[0], rs[1], r.getX());
      rs[1].setX (F2elm.mult (rs[1].getX(), r.getZ()));
    }

    return r;
  }
    
  
  public F2elm jInv () {
    // Computes the j-invariant of a Montgomery curve

    F2elm t0, t1, jinv;

    jinv = F2elm.sqr (a);
    t1 = F2elm.sqr (c);
    t0 = F2elm.leftShift (t1, 1);
    t0 = F2elm.sub (jinv, t0);
    t0.f2SubInPlace (t1);
    jinv = F2elm.sub (t0, t1);
    t1.f2SqrInPlace ();
    jinv.f2MultInPlace (t1);
    t0.f2LeftShiftInPlace (2);
    t1 = F2elm.sqr (t0);
    t0.f2MultInPlace (t1);
    t0 = F2elm.leftShift (t0, 2);
    jinv.f2InverseInPlace ();
    jinv.f2MultInPlace (t0);

    return jinv;
  }
    

  public String toString() {
    F2elm adivc;
    
    adivc = F2elm.inverse (c);
    adivc.f2MultInPlace (a);
    
    return "y^2 = x^3 + " + adivc + " x^2 + " + c + " x";
  }
}
