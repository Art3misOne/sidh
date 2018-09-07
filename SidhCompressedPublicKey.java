
package sidh;

/**************************************************************************************************
 *
 * Implements compressed public keys for the Supersingular Isogeny Diffie-Hellman key exchange 
 * algorithm. 
 *  
 **************************************************************************************************

import java.math.BigInteger;
import java.util.Arrays;

class SidhCompressedPublicKey {
  F2elm a3, a4;
  BigInteger[] triple3;
  BigInteger[] triple4;
  boolean swapcd3, swapcd4;

    
  public SidhCompressedPublicKey (SidhPublicKey pk, BigInteger order3, BigInteger order4) {
    compress3 (pk.get4Key (), pk.getP4 (), pk.getQ4 (), order4);
    compress4 (pk.get3Key (), pk.getP3 (), pk.getQ3 (), order3);
  }


  private void compress3 (F2elm a, F2elm phiPx, F2elm phiQx, BigInteger order4) {
    F2Point phiP, phiQ, s1, s2;
    F2Point[] ss;
    F2elm zInv;
    BigInteger[] cd;
    
    phiP = new F2Point (phiPx, F2elm.ONE);
    phiP.recoverYcoord (a);
    phiQ = new F2Point (phiQx, F2elm.ONE);
    phiQ.recoverYcoord (a);

    ss = get2TorsionEntangledBasis (a);
    s1 = ss[0].normalize();
    s2 = ss[1].normalize();
 
    cd = poligHelman (a, phiP, phiQ, s1, s2);
    if (cd[3].testBit (0)) {
      triple3 = computeTriple (cd[0], cd[1], cd[2], cd[3], order4);
      swapcd3 = false;
    }
    else {
      triple3 = computeTriple (cd[2], cd[3], cd[0], cd[1], order4);
      swapcd3 = true;
    }
  }


  private void compress4 (F2elm a, F2elm phiPx, F2elm phiQx, BigInteger order3) {
    F2Point phiP, phiQ, p, q;
    F2Point[] pq;
    F2elm zInv;
    BigInteger d1mod3, three = BigInteger.valueOf (3);
    BigInteger[] cd;
    
    phiP = new F2Point (phiPx, F2elm.ONE);
    phiP.recoverYcoord (a);
    phiQ = new F2Point (phiQx, F2elm.ONE);
    phiQ.recoverYcoord (a);

    pq = buildOrdinaryE3nBasis (a);
    p = pq[0].normalize ();
    q = pq[1].normalize ();

    cd = poligHelman (a, phiP, phiQ, p, q);

    d1mod3 = cd[3].mod (three);
    if (d1mod3.equals (BigInteger.ZERO)) {
      triple4 = computeTriple (cd[2], cd[3], cd[0], cd[1], order3);
      swapcd4 = false;
    }
    else {
      triple4 = computeTriple (cd[0], cd[1], cd[2], cd[3], order3);
      swapcd4 = true;
    }
  }


  private BigInteger[] computeTriple (BigInteger c0, BigInteger c1,
				      BigInteger d0, BigInteger d1, BigInteger order) {
    BigInteger c0inv;
    BigInteger[] result = new BigInteger[3];
    
    c0inv = c0.modInverse (order);
    result[0] = d1.multiply (c0inv);
    result[0] = result[0].mod (order);

    result[1] = c1.multiply (c0inv);
    result[1] = result[0].mod (order);

    result[2] = d0.multiply (c0inv);
    result[2] = result[0].mod (order);

    return result;
  }
    

  public F2Point[] decompress (BigInteger order3, BigInteger order4, SidhPrivateKey sk) {
    F2Point[] r3r4;

    r3r4[0] = decompress3 (sk, order4);
    r3r4[1] = decompress4 (sk, order3);

    return r3r4;
  }


  private F2Point decompress3 (SidhPrivateKey seckey, BigInteger order4) {
    F2Point[] ss;
    F2Point s1, s2;
    Felm a24;
    BigInteger t1, t2, t3, m4;
	
    ss = get2TorsionEntangledBasis ();
    s1 = ss[0].normalize ();
    s2 = ss[1].normalize ();
    
    // a24 = (a+2) / 4
    a24 = a4.f2Add (Felm.ONE);
    a24 = a24.f2Add (Felm.ONE);
    a24 = a24.f2Div2 ();
    a24 = a24.f2Div2 ();

    m4 = secket.getKey4 ();
    
    if (swapcd3) {
      t1 = m4.multiply (triple3[1]);
      t1 = t1.add (BigInteger.ONE);
      t2 = t1.modInverse (order4);
      t1 = m4.multiply (triple3[2]);
      t1 = t1.add (triple3[0]);
      t3 = t1.multiply (t2);
      t3 = t3.and (mask);
      r = twoDimScalarMult (t3, s1, s2, a, a24);
    }
    else {
      t1 = m4.multiply (triple3[2]);
      t1 = t1.add (BigInteger.ONE);
      t1 = t1.and (mask);
      t2 = t1.modInverse (order4);
      t1 = m4.multiply (triple3[1]);
      t1 = t1.add (BigInteger.ONE);
      t3 = t1.multiple (t2);
      t3 = t3.and (mask);
      r = twoDimScalarMult (t3, s1, s2, a, a24);
    }

    //...???
    return r;
  }


  private F2Point decompress4 (SidhPrivateKey seckey, BigInteger order3) {
    F2Point[] rs;
    F2Point r, r1, r2;
    Felm a24;
    BigInteger t1, t2, t3;
    
    rs = buildOrdinaryE3nBasis (a);
    r1 = pq[0].normalize ();
    r2 = pq[1].normalize ();

    // a24 = (a+2) / 4
    a24 = a4.f2Add (Felm.ONE);
    a24 = a24.f2Add (Felm.ONE);
    a24 = a24.f2Div2 ();
    a24 = a24.f2Div2 ();

    t1 = seckey.get3Key ();
    t2 = triple4[0];
    t3 = triple4[1];
    t4 = triple4[2];

    if (swapcd4 == false) {
      t3 = t1.multiply (t3);
      t3 = t3.add (BigInteger.ONE);
      t3 = t3.modInverse (order3);
      t4 = t1.multiply (t4);
      t4 = t2.add (t4);
      t3 = t3.multiply (t4);
      t3 = t3.mod (order3)
      r = twoDimScalarMult (t3, r1, r2, a4, a24);
    }
    else {
      t4 = t1.multiply (t4);
      t4 = t4.add (BigInteger.ONE);
      t4 = t4.modInverse (order3);
      t3 = t1.multiply (t3);
      t3 = t2.add (t3);
      t3 = t3.multiply (t4);
      t3 = t3.mod (order3);
      r = twoDimScalarMult (t3, r1, r2, a4, a24);
    }

    return r;
  }


  private F2Point[] get2TorsionEntangledBasis (int eA, int eB) {
    int i, index = 0, index2 = 0;
    Felm r, s, z, apha, twoalphainv, beta;
    F2elm t, u, u0;
    F2elm[] table;
    boolean isSqrA = false;
    F2Point[] result;
    
    if (a3.isQuadraticResidue ()) {
      table = tableVqnr;
      isSqrA = true;
    }
    else
      table = tableVqr;

    do {
      x1 = a3.f2Mult (table[index2]);
      index2++;
      x1 = x1.f2Negate ();
      t = x1.f2Add (a3);
      t = x1.f2Mult (t);
      t = t.f2Add (F2elm.ONE);
      t = t.f2Mult (x1);
      index += 2;
    } while (t.isQuadraticResidue () == false);

    if (isSqrA)
      tableVqnr[(index-2)/2] = r;
    else
      tableVqr[(index-2)/2)] = r;

    z = s.fpAdd (t.f2Get0 ());
    z = z.fpDiv2 ();
    alpha = new Felm (z);
    for (i = 0; i < (eA - 2); i++) 
      alpha = alpha.fpSqr ();
    for (i = 0; i < eB; i++) {
      s = alpha.fpSqr ();
      alpha = alpha.fpMult (s);
    }

    twoalphainv = alpha.fpAdd (alpha);
    twoalphainv = twoalphainv.fpInverse ();
    beta = twoalphainv.fpMult (t.f2Get1 ());
    twoalphainv = alpha.fpSqr ();

    if (twoalphainv.fpCompare (z)) 
      y1 = new F2elm (alpha, beta);
    else
      y1 = new F2elm (beta.fpNegate(), alpha.fpNegate());

    x2 = x1.f2Add (a3);
    x2 = x2.f2Negate ();
    y2 = u0.f2Mult (y1);
    y2 = y2.f2Mult (new F2elm (r, Felm.ZERO));

    result = new F2Point[2];
    result[0] = new F2Point (x1, y1, F2elm.ONE);
    result[1] = new F2Point (x2, y2, F2elm.ONE);

    return result;
  }
    

  /***** From Microsoft code base *****
  Felm getPointNotIn2E (Felm alpha, Felm fortyseven, Felm fiftytwo, Felm four) {
    Felm a0, a1, x0, x1, y0, y1, t0, temp0, temp1, sqrt;
    Felm alpha47, alphasq47, alpha52, alphasq52, alphaout;
    int i;
    
    a0 = a.f2Get0 ();
    a1 = a.f2Get1 ();

    x0 = a0.fpSub (a1);                             // x0 = a0 - a1
    x0 = x0.fpAdd (a0);                             // x0 = x0 + a0 = 2*a0 - a1
    x0 = x0.fpLeftShift (3);                        // x0 = 8 * (2*a0 - a1)
    y0 = x0.fpSub (a0);                             // y0 = x0 - a0 = 15*a0 - 8*a1
    x1 = a0.fpAdd (a1);                             // x1 = a0 + a1
    x1 = x1.fpAdd (a1);                             // x1 = x1 + a1 = a0 + 2a1
    y1 = x1.fpLeftShift (3);                        // y1 = 8*x1 = 8*a0 + 16*a1

    alphaout = new Felm (alpha);
    alpha52 = alphaout.fpMult (fiftytwo);           // alpha52 = 52*alpha
    alphasq52 = alpha52.fpMult (alphaout);          // alphasq52 = 52*alpha^2
    alpha47 = alphaout.fpMult (fortyseven);         // alpha47 = 47*alpha
    alphasq47 = alpha47.fpMult (alphaout);          // alphasq47 = 47*alpha^2
    temp0 = y0.fpMult (alphaout);                   // temp0 = y0*alpha = (15*a0 - 8*a1) * alpha
    temp1 = y1.fpMult (alphaout);                   // temp1 = y1*alpha = (8*a0 + 16*a1) * alpha

    do {
      alphaout = alphaout.fpAdd (Felm.ONE);         // alpha += 1
      temp0 = temp0.fpMult (y0);                    // temp0 *= y0
      t0 = alpha52.fpAdd (fiftytwo);                // t0 = 52*alpha + 52
      alpha52 = t0.fpAdd (alpha52);                 // alpha52 = 2*52*alpha + 52
      alphasq52 = alphasq52.fpAdd (alpha52);        // alphasq52 = 52*alpha^2 + 2*52*alpha + 52
      alpha52 = new Felm (t0);                      // alpha52 = 52*alpha + 52
      x0 = alphasq52.fpAdd (four);                  // x0 = alphasq52 + 4
      x0 = x0.fpAdd (temp0);                        // x0 += temp0
      temp1 = temp1.fpAdd (y1);                     // temp1 += y1
      t0 = alpha47.fpAdd (fortyseven);              // t0 = 47*alpha + 47
      alpha47 = alpha47.fpAdd (t0);                 // alpha47 = 2*47*alpha + 47
      alphasq47 = alphasq47.fpAdd (alpha47);        // alphasq47 = 47*alpha^2 + 2*47*alpha + 47
      alpha47 = new Felm (t0);                      // alpha47 = 47*alpha + 47
      x1 = alphasq47.fpAdd (Felm.ONE);              // x1 = alphasq47 + 1
      x1 = x1.fpAdd (temp1);                        // x1 += temp1
      x0 = x0.fpSqr ();
      x1 = x1.fpSqr ();
      t0 = alphaout.fpSqr ();
      x0 = x0.fpAdd (x1);                           // x0 = x0^2 + x1^2
      t0 = t0.fpMult (x0);                          // t0 = alpha^2 * (x0^2 + x1^2)
      sqrt = new Felm (t0);                     
      for (i = 0; i < 371; i++)                     // sqrt = t0 ^ ((p+1) div 2)
	sqrt = sqrt.fpSqr ();
      for (i = 0; i < 239; i++) {
	x0 = sqrt.fpSqr ();
	sqrt = sqrt.fpMult (x0);
      }
    } while (sqrt.fpEquals (t0) == false);

    return alphaout;
  }


  F2Point[] generate2TorsionBasis () {
    F2Point p, q, p1, p2, r1, r2;
    Felm alpha, four, fortyseven, fiftytwo;
    F2elm x1, x2, z1, z2, px, pz, qx, qz, y1, y2, t0, t1;

    F2Point[] retval = new F2Point[2];
    
    MontCurve curve = new MontCurve (a, F2elm.ONE);
    
    four = new Felm (4);
    fortyseven = new Felm (47);
    fiftytwo = new Felm (52);
    alpha = Felm.ZERO;
    
    alpha = getPointNotIn2E (alpha, fortyseven, fiftytwo, four);
    x1 = new F2elm (alpha.fpMult (four), alpha);    // x1 = 4*alpha + alpha*i
    
    p1 = new F2Point (x1, F2elm.ONE);               // p1 = (x1, 1)
    p1 = curve.xTple (p1, 239);                     // p1 = (x1,1)*(3^239)
    p = curve.xDble (p1, 371);                      // p = (x1,1)*(3^239)*(2^371)
    
    do {
      alpha = getPointNotIn2E (alpha, fortyseven, fiftytwo, four);
      x2 = new F2elm (alpha.fpMult (four), alpha);  // x2 = 4*alpha + alpha*i
      p2 = new F2Point (x2, F2elm.ONE);             // p2 = (x2, 1)
      p2 = curve.xTple (p2, 239);                   // p2 = (x2,1)*(3^239)
      q = curve.xDble (p2, 371);                    // q = (x2,1)*(3^239)*(2^371)
      t0 = p.getX().f2Mult (q.getZ());              // t0 = px*qz
      t1 = q.getX().f2Mult (p.getZ());              // t1 = qx*pz
      t0 = t0.f2Sub (t1);                           // t0 = px*qz - qx*pz 
    } while (t0.f2Equals (F2elm.ZERO));
    
    r1 = new F2Point (x1, z1);
    r2 = new F2Point (x2, z2);

    r1.recoverYcoord ();
    r2.recoverYcoord ();

    retval[0] = r1;
    retval[1] = r2;

    return retval;
  }


  F2elm getXonCurve (int r, Felm t1, Felm aprime, Felm b) {
    Felm v0, v1, r0, r1, t0, t2, t3, rsq, a0, a1;
    F2elm x;
    int i;

    a0 = a.f2Get0 ();
    a1 = a.f2Get1 ();
    
    // Figure out what LIST is, set r1 = list[2*r-1], r0 = list[2*r]

    rsq = r.fpSqr ();                               // rsq = r^2
    t0 = r1.fpMult (a1);                            // t0 = r1 * a1
    v0 = r0.fpMult (a0);                            // t0 = r0 * a0
    v0 = v0.fpSub (t0);                             // v0 = v0 - t0
    t0 = r0.fpMult (a1);                            // t0 = r0 * a1
    v1 = r1.fpMult (a0);                            // v1 = r1 * a0
    v1 = v1.fpAdd (t0);                             // v1 = v1 + t0
    t0 = v0.fpAdd (a0);                             // t0 = v0 + a0
    t1 = v1.fpAdd (a1);                             // t1 = v1 + a1
    t2 = v0.fpMult (v1);                            // t2 = v0 * v1
    t2 = t2.fpShiftLeft (1);                        // t2 = t2*2


    // ...

    return x;
  }


  FpPointAffine getPtOnCurve (int *r) {
    Felm t0, t1, t2, t3, aprime, b;
    F2elm x, y;
    
    getXonCurve (r, x, t1, aprime, b);

    t0 = a.fpAdd (t1);                              // t0 = a+t1
    t0 = t0.fpShiftRight (1);                       // t0 = t0/2
    t1 = new Felm (t0);                             // t1 = t0
    // t1 = t0^((p-3)/4)
    t3 = t0.fpMult (t1);                            // t3 = t0*t1
    t2 = t3.fpSqr ();                               // t2 = t3^2
    t1 = t1.fpShiftRight (1);                       // t1 = t1/2
    t1 = t1.fpMult (b);                             // t1 = t1*b

    if (t0.fpEquals (t2)) {
      y.f2Set0 (t3);                                // y0 = t3
      y.f2Set1 (t1);                                // y1 = t1
    }
    else {
      t3 = t3.fpNeg ();
      y.f2Set0 (t1);                                // y0 = t1
      y.f2Set1 (t3);                                // y1 = -t3
    }

    return new F2Point (x, y, F2elm.ONE);
  }


  get3TorsionElement (unsigned int* r, F2Point p, F2Point p3) {
    F2Point pp;
    F2elm a24, c24, x, z;
    Felm t0, t1, t2;

    MontCurve curve = new MontCurve (a, F2elm.ONE);
    
    getXonCurve (r, x, t0, t1, t2);
    pp = new F2Point (x, F2elm.ONE);
    pp = curve.xDble (pp, 372);

    pp = new F2Point (p);

    c24 = new Felm (BigInteger.valueOf (2));        // c24 = 2
    a24 = a.f2Add (c24);                            // a24 = a + 2
    c24 = c24.f2LeftShift (1);                      // c24 = 4

    z = pp.getZ ();
    while (z.f2Equals (F2elm.ZERO) == false) {      // while pp.Z != 0
      p3 = new F2Point (pp);
      pp = curve.xTpl (pp);
      triples++;
      z = pp.getZ ();
    }

    // return p3 and triples?
  }


  generate3TorsionBasis (int eA, int eB) {
    F2Point r, r1, r2, r3, r4;
    MontCurve curve = new MontCurve (a, F2elm.ONE);
    F2elm u, v, c, f, t0, f0, fx, fy, y3, x1, z1;
    int triples = 0, ptsFound = 0, ri = 1;;

    // arguments and return value?
    get3TorsionElement ();

    if (triples == eB) {
      ptsFound = 1;
      x1 = r.getX ();
      z1 = r.getZ ();
      u = a.f2Mult (z1);
      u = u.f2Add (x1);
      u = u.f2Mult (x1);
      v = z1.f2Sqr ();
      u = u.f2Add (v);
      u = u.f2Mult (x1);
      v = v.f2Mult (z1);
      y1 = u.f2Mult (v.f2Inverse());
      y1 = y1.f2Sqrt ();
      y1 = y1.f2Mult (z1);
    }

    x3 = r3.getX ();
    z3 = r3.getZ ();
    
    u = a.f2Mult (z3);
    u = u.f2Add (x3);
    u = u.f2Mult (x3);
    v = z3.f2Sqr ();
    u = u.f2Add (v);
    u = u.f2Mult (x3);
    v = v.f2Mult (z3);
    y3 = u.f2Mult (v.f2Inverse ());
    y3 = y3.f2Sqrt ();
    y3 = y3.f2Mult (z3);
    f0 = x3.f2Sqr ();
    t0 = z3.f2Sqr ();
    fx = x3.f2Mult (z3);
    fx = a.f2Mult (fx);
    fx = fx.f2LeftShift (1);
    fx = fx.f2Add (t0);
    fx = fx.f2Add (f0);
    fx = fx.f2Add (f0);
    fx = fx.f2Add (f0);
    f0 = t0.f2Sub (f0);
    fx = fx.f2Mult (z3);
    fy = y3.f2Mult (z3);
    fy = fy.f2LeftShift (1);
    fy = fy.f2Negate ();
    c = fy.f2LeftShift (1);
    fy = fy.f2Mult (z3);
    f0 = f0.f2Mult (x3);
    c = c.f2Mult (y3);
    fx = c.f2Mult (fx);
    fy = c.f2Mult (fy);
    f0 = c.f2Mult (f0);

    do {
      while (ptsFound < 2) {
	ri++;
	getPtOnCurve ();
	
      }
    }
  }


  dblAndLine (F2Point p, F2el lx, F2elm ly, F2elm l0, f2elm v0) {
    F2elm x2, xz, yz, z2;
    F2elm xx2, t0;

    xx2 = yz.f2LeftShift (1);                       // xx2 = 2 * yz
    ly = xx2.f2Sqr ();                              // ly = xx2^2
    l0 = x2.f2Sub (z2);                             // l0 = x2 - z2
    v0 = l0.f2Sqr ();                               // v0 = l0^2
    l0 = xx2.f2Mult (l0);                           // l0 = xx2 * l0
    lx = xz.f2Mult (l0);                            // lx = xz * l0
    xx2 = yz.f2Mult (ly);                           // xx2 = yz * ly
    lx = xx2.f2Add (lx);                            // lx = xx2 + lx
    yz = x2.f2Add (z2);                             // yz = x2 + z2
    yz = a.f2Mult (yz);                             // yz = a * yz
    xx2 = xz.f2LeftShift (1);                       // xx2 = 2 * xz
    yz = xx2.f2Add (yz);                            // yz = xx2 + yz
    yz = xx2.f2Add (yz);                            // yz = xx2 + yz
    yz = xx2.f2Mult (yz);                           // yz = xx2 * yz

    xx2 = v0.f2Sqr ();                              // xx2 = v0^2
    t0 = l0.f2Sqr ();                               // t0 = l0^2
    z2 = ly.f2Sqr ();                               // z2 = ly^2
    yz = v0.f2Add (yz);                             // yz = v0 + yz
    yz = l0.f2Mult (yz);                            // yz = l0 * yz

    ly = xz.f2Mult (ly);                            // ly = xz * ly
    l0 = x2.f2Mult (l0);                            // l0 = x2 * l0
    v0 = xz.f2Mult (v0);                            // v0 = xz * v0

    x2 = new F2elm (xx2);
    xz = new F2elm (xz);
  }


  absorbLine (F2elm lx, F2elm ly, F2elm l0, F2elm v0, F2Point p, F2elm n, F2elm d) {
    F2elm x, y, l, v;

    x = p.getX ();
    y = p.getY ();
 
    l = lx.f2Mult (x);                              // l = lx * x
    v = ly.f2Mult (y);                              // v = ly * y
    l = v.f2Sub (l);                                // l = v - l
    l = l.f2Add (l0);                               // l = l + l0
    v = ly.f2Mult (x);                              // v = ly * x
    v = v.f2Add (v0);                               // v = v + v0
    n = n.f2Mult (l);                               // n = n * l
    d = d.f2Mult (v);                               // d = d * v

    // return (n, d)?
  }


  squareAndAbsorbLine (F2elm lx, F2elm ly, F2elm l0, F2elm v0, F2Point p, F2elm n, F2elm d) {
    n = n.f2Sqr ();
    d = d.f2Sqr ();
    absorb_line (lx, ly, l0, v0, p, n, d);

    // return (n,d)?
  }


  finalDblIteration (F2Point p, F2elm x, F2elm n, F2elm d) {
    F2elm px, pz, l;

    px = p.getX ();
    pz = p.getZ ();

    n = n.f2Sqr ();                                 // n = n^2
    d = d.f2Sqr ();                                 // d = d^2
    d = d.f2Mult (pz);                              // d = d * pz
    l = pz.f2Mult (x);                              // l = pz * x
    l = l.f2Sub (px);                               // l = l - px
    n = n.f2Mult (l);                               // n = n * l

    // return (n, d)?
  }


  F2elm finalExponentiation2Torsion (F2elm n, F2elm d, F2elm ninv, F2elm dinv) {
    int i;
    MontCurve curve = new MontCurve (a);

    n = n.f2Mult (dinv);                            // n = n * dinv
    // n = n^p ... ???
    d = d.f2Mult (ninv);                            // d = d * ninv
    n = n.f2Mult (d);                               // n = n * d

    for (i = 0; i < 239; i++) {
      // ... ???
    }

    return n;
  }


  tatePairings2Torsion (F2Point R1, F2Point R2, F2point p, F2Point q, F2elm n) {
    
    MontCurve curve = new MontCurve (a);
  }
}
*/

