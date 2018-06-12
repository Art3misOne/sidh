
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
    coeff = new F2elm[5];

    a = ia;
    b = F2elm.ONE;
    c = F2elm.ONE;

    updateA24C24();
  }


  public FourIsogeny (MontCurve curve) {
    coeff = new F2elm[5];
    
    a = curve.a;
    b = curve.b;
    c = curve.c;
    a24 = curve.a24;
    c24 = curve.c24;
  }


  public F2elm[] getCoeff() {
    return coeff;
  }


  public F2Point first4Isog (F2Point p) {
    // Compute first 4-isogeny computed by party A.
    // Inputs: projective point (px:pz)
    // Outputs: updates a and c to isogenous curve parameters and returns a point in the 
    //          codomain 

    F2elm t0, t1, t2, px, pz;
    int i;

    px = p.getX();
    pz = p.getZ();

    t0 = F2elm.ONE;
    t0 = t0.f2Add(t0);
    c = a.f2Sub(t0);                     // c = a - 2

    t0 = new F2elm (6, 0);    
    a = a.f2Add(t0);                     // a = a + 6
    a = a.f2Add(a);                      // a = 2*a + 12
    t1 = px.f2Add(pz);                   // t1 = px + pz
    t2 = px.f2Sub(pz);                   // t2 = px - pz
    t1 = t1.f2Sqr();                     // t1 = (px+pz)^2
    t2 = t2.f2Sqr();                     // t2 = (px-pz)^2
    pz = px.f2Mult(pz);                  // pz = px * pz
    pz = pz.f2Negate();                  // pz = -1 * px * pz
    pz = pz.f2Mult(c);                   // pz = -1*px*pz*c
    px = t1.f2Sub(pz);                   // px = (px+pz)^2 + c*px*pz
    pz = t2.f2Mult(pz);                  // pz = (px-pz)^2 * (-c*px*pz)
    px = px.f2Mult(t1);                  // px = (px+pz)^2 * [(px+pz)^2+c*px*pz]

    updateA24C24();

    return new F2Point (px, pz);
  }


  public void get4Isog (F2Point p) {
    // Computes the corresponding 4-isogeny of a projective Montgomery point (px:pz) of order 4.
    // Inputs: Projective point (px:pz) of order 4  
    // Outputs: updates curve parameters (a and c) and coefficients to evaluate the isogeny 

    F2elm px, pz;

    px = p.getX();
    pz = p.getZ();

    coeff[0] = px.f2Add(pz);             // coeff[0] = px + pz
    coeff[3] = px.f2Sqr();               // coeff[3] = px^2
    coeff[4] = pz.f2Sqr();               // coeff[4] = pz^2
    coeff[0] = coeff[0].f2Sqr();         // coeff[0] = (px+pz)^2
    coeff[1] = coeff[3].f2Add(coeff[4]); // coeff[1] = px^2 + pz^2
    coeff[2] = coeff[3].f2Sub(coeff[4]); // coeff[2] = px^2 - pz^2
    coeff[3] = coeff[3].f2Sqr();         // coeff[3] = px^4
    coeff[4] = coeff[4].f2Sqr();         // coeff[4] = pz^4
    coeff[0] = coeff[0].f2Sub(coeff[1]); // coeff[0] = (px+pz)^2 - (px^2 + pz^2) = 2*px*pz
    a = coeff[3].f2Add(coeff[3]);        // a = 2*px^4
    a = a.f2Sub(coeff[4]);               // a = 2*px^4 - pz^4
    a = a.f2Add(a);                      // a = 2*(2*px^4 - pz^4)
    c = coeff[4];                        // c = pz^4

    updateA24C24();
  }


  public F2Point eval4Isog (F2Point p) {
    // Evaluate the isogeny at the point (px:pz) 
    // Inputs: projective point (px:pz)
    // Outputs: Returns the corresponding point in the domain of the isogeny

    F2elm t0, t1, px, pz;

    px = p.getX();
    pz = p.getZ();

    px = px.f2Mult(coeff[0]);            // px = px * coeff[0]
    t0 = pz.f2Mult(coeff[1]);            // t0 = pz * coeff[1]
    px = px.f2Sub(t0);                   // px = coeff[0]*px - coeff[1]*pz
    pz = pz.f2Mult(coeff[2]);            // pz = coeff[2] * pz
    t0 = px.f2Sub(pz);                   // t0 = coeff[0]*px - coeff[1]*pz - coeff[2]*pz
    pz = pz.f2Mult(px);                  // pz -> px * pz
    t0 = t0.f2Sqr();                     // t0 -> t0^2
    pz = pz.f2Add(pz);                   // pz -> pz+pz
    pz = pz.f2Add(pz);                   // pz -> pz+pz
    px = pz.f2Add(t0);                   // px -> t0+pz
    pz = pz.f2Mult(t0);                  // pz -> t0*pz
    pz = pz.f2Mult(coeff[4]);            // pz -> pz*coeff[4]
    t0 = t0.f2Mult(coeff[4]);            // t0 -> t0*coeff[4]
    t1 = px.f2Mult(coeff[3]);            // t1 -> px * coeff[3]
    t0 = t0.f2Sub(t1);                   // t0 -> t0 - t1
    px = px.f2Mult(t0);                  // px -> px*t0

    return new F2Point (px, pz);
  }
}


class ThreeIsogeny extends MontCurve {

  public ThreeIsogeny (F2elm ia) {
    a = ia;
    b = F2elm.ONE;
    c = F2elm.ONE;

    updateA24C24();
  }


  public ThreeIsogeny (MontCurve curve) {
    a = curve.a;
    b = curve.b;
    c = curve.c;
    a24 = curve.a24;
    c24 = curve.c24;
  }


  public void get3Isog (F2Point p) {
    // Computes the corresponding 3-isogeny of a projective Montgomery point of order 3 (px:pz)
    // Inputs: Projective point of order 3 (px:pz)
    // Outputs: Computes curve parameters a, c 

    F2elm t0, t1, px, pz;

    px = p.getX();
    pz = p.getZ();

    t0 = px.f2Sqr();                     // t0 = x^2
    t1 = t0.f2Add(t0);                   // t1 = 2*t0
    t0 = t0.f2Add(t1);                   // t0 = t0+t1
    t1 = pz.f2Sqr();                     // t1 = z^2
    a = t1.f2Sqr();                      // a = t1^2
    t1 = t1.f2Add(t1);                   // t1 = 2*t1
    c = t1.f2Add(t1);                    // c = 2*t1
    t1 = t0.f2Sub(t1);                   // t1 = t0-t1
    t1 = t0.f2Mult(t1);                  // t1 = t0*t1
    a = a.f2Sub(t1);                     // a = a-t1 
    a = a.f2Sub(t1);                     // a = a-t1 
    a = a.f2Sub(t1);                     // a = a-t1     
    t1 = px.f2Mult(pz);                  // t1 = x*z    
    c = c.f2Mult(t1);                    // c = c*t1    

    updateA24C24();
  }


  public F2Point eval3Isog (F2Point p, F2Point q) {
    // Evaluate the isogeny at q = (x:z)
    // Inputs: projective point p = (x3:z3) of order 3 and projective point q
    // Outputs: Returns phi(q)

    F2elm t0, t1, t2, px, pz, qx, qz;

    px = p.getX();
    pz = p.getZ();

    qx = q.getX();
    qz = q.getZ();

    t0 = px.f2Mult(qx);                  // t0 = x3*x
    t1 = pz.f2Mult(qx);                  // t1 = z3*x
    t2 = pz.f2Mult(qz);                  // t2 = z3*z
    t0 = t0.f2Sub(t2);                   // t0 = x3*x - z3*z
    t2 = px.f2Mult(qz);                  // t2 = x3*z
    t1 = t1.f2Sub(t2);                   // t1 = z3*x - x3*z
    t0 = t0.f2Sqr();                     // t0 = (x3*x - z3*z)^2
    t1 = t1.f2Sqr();                     // t1 = (z3*x - x3*z)^2
    qx = qx.f2Mult(t0);                  // x = x*(x3*x - z3*z)^2
    qz = qz.f2Mult(t1);                  // z = z*(z3*x - x3*z)^2

    return new F2Point (qx, qz);
  }
}


