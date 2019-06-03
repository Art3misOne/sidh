
package sidh;

/**************************************************************************************************
 *
 * Implements classes for elliptic curve isogenies of degree three and four. These extend the 
 * MontCurve class as an isogeny is a curve with additional methods for computing and evaluating
 * the isogeny.
 *  
 **************************************************************************************************/


import java.math.BigInteger;


class FourIsogeny extends MontCurve {
  F2elm coeff[];


  public FourIsogeny (F2elm ia) {
    coeff = new F2elm[3];

    a = ia;
    c = F2elm.ONE;
  }


  public FourIsogeny (MontCurve curve) {
    coeff = new F2elm[3];
    
    a = curve.a;
    c = curve.c;
    a24 = curve.a24;
    c4 = curve.c4;
    aPlus2c = curve.aPlus2c;
    aMinus2c = curve.aMinus2c;
  }


  public F2elm[] getCoeff() {
    return coeff;
  }


  public void get4Isog (F2Point p) {
    // Compute 4-isogeny of a point (px:pz) of order 4. Updates curve parameters and coefficients. 

    F2elm px, pz;

    px = p.getX();
    pz = p.getZ();

    coeff[1] = F2elm.sub (px, pz);
    coeff[2] = F2elm.add (px, pz);
    coeff[0] = F2elm.sqr (pz);
    coeff[0].f2LeftShiftInPlace (1);
    c4 = F2elm.sqr (coeff[0]);
    coeff[0].f2LeftShiftInPlace (1);
    aPlus2c = F2elm.sqr (px);
    aPlus2c.f2LeftShiftInPlace (1);
    aPlus2c.f2SqrInPlace ();
  }


  public F2Point eval4Isog (F2Point p) {
    F2elm t0, t1, px, pz;

    px = p.getX ();
    pz = p.getZ ();

    t0 = F2elm.add (px, pz);
    t1 = F2elm.sub (px, pz);
    px = F2elm.mult (t0, coeff[1]);
    pz = F2elm.mult (t1, coeff[2]);
    t0.f2MultInPlace (t1);
    t0.f2MultInPlace (coeff[0]);
    t1 = F2elm.add (px, pz);
    pz = F2elm.sub (px, pz);
    t1.f2SqrInPlace ();
    pz.f2SqrInPlace ();
    px = F2elm.add (t0, t1);
    t0 = F2elm.sub (pz, t0);
    px.f2MultInPlace (t1);
    pz.f2MultInPlace (t0);
    
    return new F2Point (px, pz);
  }
}


class ThreeIsogeny extends MontCurve {
  F2elm coeff[];

    
  public ThreeIsogeny (F2elm ia) {
    coeff = new F2elm[2];

    a = ia;
    c = F2elm.ONE;
  }


  public ThreeIsogeny (MontCurve curve) {
    coeff = new F2elm[2];
      
    a = curve.a;
    c = curve.c;
    a24 = curve.a24;
    c4 = curve.c4;
    aPlus2c = curve.aPlus2c;
    aMinus2c = curve.aMinus2c;
  }


  public F2elm[] getCoeff() {
    return coeff;
  }


  public void get3Isog (F2Point p) {
    // Compute 3-isogeny of a point (px:pz) of order 3. Updates curve parameters and coefficients.

    F2elm t0, t1, t2, t3, t4, px, pz;

    px = p.getX ();
    pz = p.getZ ();

    coeff[0] = F2elm.sub (px, pz);
    t0 = F2elm.sqr (coeff[0]);
    coeff[1] = F2elm.add (px, pz);
    t1 = F2elm.sqr (coeff[1]);
    t2 = F2elm.add (t0, t1);
    t3 = F2elm.add (coeff[0], coeff[1]);
    t3.f2SqrInPlace ();
    t3.f2SubInPlace (t2);
    t2 = F2elm.add (t1, t3);
    t3.f2AddInPlace (t0);
    t4 = F2elm.add (t0, t3);
    t4.f2LeftShiftInPlace (1);
    t4.f2AddInPlace (t1);
    aMinus2c = F2elm.mult (t2, t4);
    t4 = F2elm.add (t1, t2);
    t4.f2LeftShiftInPlace (1);
    t4.f2AddInPlace (t0);
    t4.f2MultInPlace (t3);
    t0 = F2elm.sub (t4, aMinus2c);
    aPlus2c = F2elm.add (t0, aMinus2c);
  }


  public F2Point eval3Isog (F2Point q) {
    // Evaluate the isogeny at q = (x:z)

    F2elm t0, t1, t2, qx, qz;

    qx = q.getX ();
    qz = q.getZ ();

    t0 = F2elm.add (qx, qz);
    t1 = F2elm.sub (qx, qz);
    t0.f2MultInPlace (coeff[0]);
    t1.f2MultInPlace (coeff[1]);
    t2 = F2elm.add (t1, t0);
    t0 = F2elm.sub (t1, t0);
    t2.f2SqrInPlace ();
    t0.f2SqrInPlace ();
    qx.f2MultInPlace (t2);
    qz.f2MultInPlace (t0);

    return new F2Point (qx, qz);
  }
}


