
package sidh;

/**************************************************************************************************
 *
 * Implements validation of public keys for the supersingular isogeny Diffie-Hellman key exchange
 * algorithm. 
 * 
 * Failure to validate public keys may result in a user generating a malicious public key that 
 * allowing them to ascertain information about a private key. However, validation of public keys
 * is more computationally expensive than generating keys. To ensure security without requiring
 * validation of public keys, private keys should only be used for a single exchange. If it is 
 * impractical for a particular application to generate a new private key for each exchange then
 * the below validation methods should be used to protect private keys.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Random;

class SidhKeyValidate {
  SidhPublicKey pk;


  public SidhKeyValidate (SidhPublicKey k) {
    pk = new SidhPublicKey (k);
  }
    

  private BigInteger genRandom (BigInteger bound) {
    // Uses BigInteger's random constructor to generate random values up to the same bitlength 
    // as the bound until the value generated is strictly less than the bound. Average expected
    // number of calls is less than 2.

    // Note: the random generator need not be cryptographically secure in this context since it is
    // only used for probabilistically checking that the curve is in the right supersingular isogeny 
    // class by verifying that a random point has the correct order.
    
    Random rnd = new Random();
    int numBits = bound.bitLength();
    
    BigInteger randval = new BigInteger (numBits, rnd);
        
    while (randval.compareTo(bound) >= 0)
      randval = new BigInteger (numBits, rnd);

    return randval;
  }


  private F2elm genRandom() {
    BigInteger x0, x1;

    x0 = genRandom (Felm.getPrime());
    x1 = genRandom (Felm.getPrime());

    return new F2elm (x0, x1);
  }


  public Boolean testCurve (int eA, int eB) {
    // Verify:
    //    1. That the the curve is not a subfield curve.
    //    2. That the curve is in the correct supersingular isogeny class.

    Felm c0d1, c1d0;
    F2Point rp, p1;
    F2elm t0, t1, rx, rz, px, pz, a = pk.getA();
    MontCurve curve = new MontCurve (a, F2elm.ONE);

    // Check the 1st condition by verifying that the j-invariant is not in Fp without 
    // computing the full j-invariant
    t0 = a.f2Sqr();                               // t0 = a^2
    t0 = t0.f2Sub (F2elm.ONE);
    t0 = t0.f2Sub (F2elm.ONE);
    t0 = t0.f2Sub (F2elm.ONE);                     // t0 = a^2 - 3
    t1 = t0.f2Sqr();                              
    t1 = t1.f2Mult (t0);                           // t1 = (a^2 - 3)^3 = t0^3
    t0 = t0.f2Sub (F2elm.ONE);                     // t0 = a^2 - 4
    c0d1 = t0.f2Get0().fpMult(t1.f2Get1());        // c0d1 = t0.x0 * t1.x1
    c1d0 = t0.f2Get1().fpMult(t1.f2Get0());        // c1d0 = t1.x0 * t0.x1
    
    if (c0d1.fpEquals (c1d0))
      return false;

    // Check the 2nd (probabilistically) by verifying that a random point has order p-1 or p+1
    rx = genRandom();
    rp = new F2Point (rx, F2elm.ONE);              // rp = (rx:1)
    rp = curve.xDble (rp, 1);                      // rp = 2*(rx:1)
    p1 = curve.xDble (rp, eA-1);                   // p1 = 2^eA*(rx:1)
    p1 = curve.xTple (p1, eB);                     // p1 = 2^eA*3^eB(rx:1)

    px = p1.getX();
    pz = p1.getZ();

    rx = pz.f2Mult (rp.getX());                    // rx = rpx * pz
    rz = px.f2Mult (rp.getZ());                    // rz = rpz * px
    rx = rx.f2Sub (rz);                            // rx = rpx*pz - rpz*px
    rx = rx.f2Mult (pz);                           // rx = pz * (rpx*pz - rpz*px)

    return rx.f2Equals (F2elm.ZERO);
  }


  public static F2elm[] cubeIndeterminant (F2elm a, F2elm b, F2elm ysq) {
    // Find na and nb such that na*y+nb = (a*y+b)^3

    F2elm t0, t1, t2, t3, ab[] = new F2elm[2];

    t0 = a; 
    t1 = b.f2Sqr();                               // t1 = b^2
    t2 = t0.f2Sqr();                               // t2 = a^2
    t2 = ysq.f2Mult (t2);                          // t2 = a^2 * ysq
    t3 = t1.f2Add (t2);                            // t3 = b^2 + a^2*ysq
    ab[0] = t1.f2Add (t3);                         // a = 2b^2 + a^2*ysq
    ab[0] = t1.f2Add (ab[0]);                      // a = 3b^2 + a^2*ysq
    ab[0] = t0.f2Mult (ab[0]);                     // a = a * (3b^2 + a^2*ysq)
    t1 = t2.f2Add (t3);                            // t1 = b^2 + 2*a^2*ysq
    t1 = t2.f2Add (t1);                            // t1 = b^2 + 3*a^2*ysq
    ab[1] = b.f2Mult (t1);                         // b = b * (b^2 + 3*a^2*ysq)

    // na*y+nb = 3*a*b^2*y + a^3*y^3 + b^3 + 3*a^2*b*y^2 = (a*y+b)^3

    return ab;
  }


  public static F2elm[] lineIndeterminantTpl (F2elm a, F2elm b, F2elm c, F2elm d, F2elm ysq) {
    // Find na and nb such that na*y+nb = (a*y+b)*(c*y+d)
    F2elm t0, t1, ab[] = new F2elm[2];

    t0 = a.f2Mult (c);                             // t0 = a*c
    ab[0] = a.f2Mult (d);                          // na = a*d
    t1 = b.f2Mult (c);                             // t1 = b*c
    ab[0] = ab[0].f2Add (t1);                      // na = a*d + b*c
    t1 = b.f2Mult (d);                             // t1 = b*d
    ab[1] = t0.f2Mult (ysq);                       // nb = a*c*y^2
    ab[1] = ab[1].f2Add (t1);                      // nb = a*c*y^2 + b*d

    // na*y+nb = a*d*y + b*c*y + a*c*y^2 + b*d = (a*y + b)*(c*y+d) 

    return ab;
  }


  public F2elm[] tplLine (F2Point p, F2Point q) {
    F2elm retvals[] = new F2elm[12];
    F2elm x3, z3, x4, z4, xUP, zUP, xUQ, zUQ, aNumer, aDenom, bNumer, bDenom;
    F2elm px, pz, qx, qz, a, xP, xQ; 
    F2elm t0, t1, t2, t3, t4, t5, t6, l0p, l1p, l2p, l0q, l1q, l2q;
    
    px = p.getX();  pz = p.getZ();
    qx = q.getX();  qz = q.getZ();

    a = pk.getA();
    xP = pk.getP();
    xQ = pk.getQ();

    t0 = px.f2Sqr();                              // t0 = px^2
    t1 = pz.f2Sqr();                              // t1 = pz^2
    t2 = px.f2Mult (pz);                           // t2 = px*pz
    t3 = a.f2Mult (t2);                            // t3 = a*t2
    t4 = t0.f2Add (t1);                            // t4 = t0+t1
    t5 = t3.f2Add (t4);                            // t5 = t4+t3
    t3 = t3.f2Add (t5);                            // t3 = t3+t5
    t5 = t5.f2Add (t5);                            // t5 = t5+t5
    l2p = t1.f2Add (t1);                           // l2p = t1+t1
    l2p = t1.f2Add (l2p);                          // l2p = l2p+t1
    l2p = t3.f2Add (l2p);                          // l2p = l2p+t3
    l2p = t5.f2Add (l2p);                          // l2p = l2p+t5
    l2p = t0.f2Mult (l2p);                         // l2p = l2p*t0
    aNumer = t1.f2Sqr();                          // aNumer = t1^2
    l2p = l2p.f2Sub (aNumer);                      // l2p = l2p-aNumer
    l1p = t0.f2Add (t3);                           // l1p = t0+t3
    l1p = t0.f2Add (l1p);                          // l1p = l1p+t0
    l1p = t5.f2Mult (l1p);                         // l1p = t5*l1p
    l1p = l1p.f2Sub (l2p);                         // l1p = l1p-l2p
    l1p = l1p.f2Add (l1p);                         // l1p = l1p+l1p
    l0p = t0.f2Sub (t1);                           // l0p = t0-t1
    l0p = t5.f2Mult (l0p);                         // l0p = l0p*t5
    l0p = l0p.f2Add (l0p);                         // l0p = l0p+l0p
    l0p = l2p.f2Sub (l0p);                         // l0p = l2p-l0p
    x3 = l0p.f2Sqr();                             // x3 = l0p^2
    z3 = l2p.f2Sqr();                             // z3 = l2p^2
    aNumer = px.f2Mult (t5);                       // aNumer = x*t5
    aNumer = aNumer.f2Add (aNumer);                // aNumer = aNumer+aNumer
    aNumer = aNumer.f2Add (aNumer);                // aNumer = aNumer+aNumer
    t0 = t0.f2Mult (l0p);                          // t0 = t0*l0p
    t5 = l2p.f2Mult (xQ);                          // t5 = l2p*xQ
    t5 = t1.f2Mult (t5);                           // t5 = t5*t1
    bNumer = l1p.f2Mult (t2);                      // bNumer = l1p*t2
    bNumer = t5.f2Add (bNumer);                    // bNumer = bNumer+t5
    bNumer = xQ.f2Mult (bNumer);                   // bNumer = bNumer*xQ
    bNumer = t0.f2Add (bNumer);                    // bNumer = bNumer+t0
    bNumer = bNumer.f2Negate();                   // bNumer = -bNumer
    t5 = a.f2Mult (t4);                            // t5 = a*t4
    t4 = t4.f2Sqr();                              // t4 = t4^2
    t2 = t2.f2Add (t2);                            // t2 = t2+t2
    xUP = t2.f2Add (t5);                           // xUP = t5+t2
    t5 = t5.f2Sub (t2);                            // t5 = t5-t2
    xUP = t2.f2Mult (xUP);                         // xUP = t2*xUP
    xUP = t4.f2Add (xUP);                          // xUP = xUP+t4
    t2 = t2.f2Mult (t5);                           // t2 = t2*t5
    t2 = t2.f2Add (t4);                            // t2 = t2+t4
    t4 = t4.f2Add (t4);                            // t4 = t4+t4
    t2 = t2.f2Add (t4);                            // t2 = t2+t4
    xUP = xUP.f2Mult (t2);                         // xUP = xUP*t2
    xUP = xUP.f2Add (xUP);                         // xUP = xUP+xUP
    xUP = xUP.f2Add (xUP);                         // xUP = xUP+xUP
    xUP = xUP.f2Sub (x3);                          // xUP = xUP-x3
    x3 = x3.f2Mult (px);                           // x3 = x*x3
    xUP = xUP.f2Sub (z3);                          // xUP = xUP-z3
    xUP = xUP.f2Negate();                         // xUP = -xUP
    xUP = l0p.f2Mult (xUP);                        // xUP = xUP*l0p
    zUP = z3.f2Mult (l2p);                         // zUP = z3*l2p
    zUP = zUP.f2Add (zUP);                         // zUP = zUP+zUP
    z3 = z3.f2Mult (pz);                           // z3 = pz*z3 
    t0 = qx.f2Sqr();                              // t0 = qx^2
    t6 = qz.f2Sqr();                              // t6 = q^2
    t2 = qx.f2Mult (qz);                           // t2 = qx*qz
    t3 = a.f2Mult (t2);                            // t3 = a*t2
    t4 = t0.f2Add (t6);                            // t4 = t0+t6
    t5 = t4.f2Add (t3);                            // t5 = t4+t3
    t3 = t3.f2Add (t5);                            // t3 = t3+t5
    t5 = t5.f2Add (t5);                            // t5 = t5+t5
    l2q = t6.f2Add (t6);                           // l2q = t6+t6
    l2q = t6.f2Add (l2q);                          // l2q = l2q+t6
    l2q = t3.f2Add (l2q);                          // l2q = l2q+t3
    l2q = t5.f2Add (l2q);                          // l2q = l2q+t5
    l2q = t0.f2Mult (l2q);                         // l2q = l2q*t0
    aDenom = t6.f2Sqr();                          // aDenom = t6^2
    l2q = l2q.f2Sub (aDenom);                      // l2q = l2q-aDenom
    l1q = t0.f2Add (t3);                           // l1q = t0+t3
    l1q = t0.f2Add (l1q);                          // l1q = l1q+t0
    l1q = t5.f2Mult (l1q);                         // l1q = t5*l1q
    l1q = l1q.f2Sub (l2q);                         // l1q = l1q-l2q
    l1q = l1q.f2Add (l1q);                         // l1q = l1q+l1q
    l0q = t0.f2Sub (t6);                           // l0q = t0-t6
    l0q = t5.f2Mult (l0q);                         // l0q = l0q*t5
    l0q = l0q.f2Add (l0q);                         // l0q = l0q+l0q
    l0q = l2q.f2Sub (l0q);                         // l0q = l2q-l0q
    x4 = l0q.f2Sqr();                             // x4 = l0q^2
    z4 = l2q.f2Sqr();                             // z4 = l2q^2
    aDenom = qx.f2Mult (t5);                       // aDenom = qx*t5
    aDenom = aDenom.f2Add (aDenom);                // aDenom = aDenom+aDenom
    aDenom = aDenom.f2Add (aDenom);                // aDenom = aDenom+aDenom
    t0 = t0.f2Mult (l0q);                          // t0 = t0*l0q
    t5 = l2q.f2Mult (xP);                          // t5 = l2q*xP
    t5 = t6.f2Mult (t5);                           // t5 = t5*t6
    bDenom = l1q.f2Mult (t2);                      // bDenom = l1q*t2
    bDenom = t5.f2Add (bDenom);                    // bDenom = bDenom+t5
    bDenom = bDenom.f2Mult (xP);                   // bDenom = bDenom*xP
    bDenom = t0.f2Add (bDenom);                    // bDenom = bDenom+t0
    bDenom = bDenom.f2Negate();                   // bDenom = -bDenom
    t5 = a.f2Mult (t4);                            // t5 = a*t4
    t4 = t4.f2Sqr();                              // t4 = t4^2
    t2 = t2.f2Add (t2);                            // t2 = t2+t2
    xUQ = t5.f2Add (t2);                           // xUQ = t5+t2
    t5 = t5.f2Sub (t2);                            // t5 = t5-t2
    xUQ = xUQ.f2Mult (t2);                         // xUQ = t2*xUQ
    xUQ = xUQ.f2Add (t4);                          // xUQ = xUQ+t4
    t2 = t2.f2Mult (t5);                           // t2 = t2*t5
    t2 = t4.f2Add (t2);                            // t2 = t2+t4
    t4 = t4.f2Add(t4);                             // t4 = t4+t4
    t2 = t2.f2Add (t4);                            // t2 = t2+t4
    xUQ = xUQ.f2Mult (t2);                         // xUQ = xUQ*t2
    xUQ = xUQ.f2Add (xUQ);                         // xUQ = xUQ+xUQ
    xUQ = xUQ.f2Add (xUQ);                         // xUQ = xUQ+xUQ
    xUQ = xUQ.f2Sub (x4);                          // xUQ = xUQ-x4
    x4 = x4.f2Mult (qx);                           // x4 = qx*x4
    xUQ = xUQ.f2Sub (z4);                          // xUQ = xUQ-z4
    xUQ = xUQ.f2Negate();                         // xUQ = -xUQ
    xUQ = l0q.f2Mult (xUQ);                        // xUQ = xUQ*l0q
    zUQ = z4.f2Mult (l2q);                         // zUQ = z4*l2q
    zUQ = zUQ.f2Add (zUQ);                         // zUQ = zUQ+zUQ
    z4 = z4.f2Mult (qz);                           // z4 = qz*z4 
    t2 = t1.f2Mult (t6);                           // t2:=t1*t6;
    t6 = t6.f2Mult (z3);                           // t6:=t6*z3;
    t1 = t1.f2Mult (z4);                           // t1:=t1*z4;
    aDenom = aDenom.f2Mult (qz);                   // aDenom:=aDenom*qz;
    aDenom = aDenom.f2Mult (z4);                   // aDenom:=aDenom*z4;
    aNumer = aNumer.f2Mult (pz);                   // aNumer:=aNumer*pz; 
    aNumer = aNumer.f2Mult (z3);                   // aNumer:=aNumer*z3;
    t3 = xP.f2Mult (z4);                           // t3:=xP*z4;
    t3 = t3.f2Sub (x4);                            // t3:=t3-x4;
    t3 = t3.f2Mult (l2q);                          // t3:=t3*l2q;
    t5 = xQ.f2Mult (z3);                           // t5:=xQ*z3;
    t5 = t5.f2Sub (x3);                            // t5:=t5-x3;
    t5 = t5.f2Mult (l2p);                          // t5:=t5*l2p;
    aNumer = aNumer.f2Mult (t3);                   // aNumer:=aNumer*t3;
    aNumer = aNumer.f2Mult (t2);                   // aNumer:=aNumer*t2;
    bNumer = bNumer.f2Mult (t3);                   // aNumer:=aNumer*t3;
    bNumer = bNumer.f2Mult (t6);                   // aNumer:=aNumer*t6;    
    aDenom = aDenom.f2Mult (t5);                   // aDenom:=aDenom*t5;
    aDenom = aDenom.f2Mult (t2);                   // aDenom:=aDenom*t2;
    bDenom = bDenom.f2Mult (t5);                   // aDenom:=aDenom*t5;
    bDenom = bDenom.f2Mult (t1);                   // aDenom:=aDenom*t2;

    retvals[0] = x3;      retvals[1] = z3;      retvals[2] = x4;       retvals[3] = z4;
    retvals[4] = xUP;     retvals[5] = zUP;     retvals[6] = xUQ;      retvals[7] = zUQ;
    retvals[8] = aNumer;  retvals[9] = aDenom;  retvals[10] = bNumer;  retvals[11] = bDenom;

    return retvals;
  } 

    
  public boolean validateAsKey (int eA, int eB) {
    F2elm t0, t1, t2, t3, t4, t5, t6, t7, lnQ, lnP, ldQ, ldP, uP, uQ, uPD, uQD, lambdaP;
    F2elm lambdaQ, sqP, sqQ, sq, alphan, betan, alphad, betad, numers[], denoms[], tlout[];
    F2elm a = pk.getA(), xP = pk.getP(), xQ = pk.getQ(), xQP = pk.getD(), px, pz, qz, qx;
    F2Point p, q, up, uq;
    int j;

    p = new F2Point (xP, F2elm.ONE);
    q = new F2Point (xQ, F2elm.ONE);

    numers = new F2elm[2];
    denoms = new F2elm[2];

    numers[0] = F2elm.ZERO;
    denoms[0] = F2elm.ZERO;
    numers[1] = F2elm.ONE;
    denoms[1] = F2elm.ONE;

    uP = F2elm.ONE;
    uQ = F2elm.ONE;
    uPD = F2elm.ONE;
    uQD = F2elm.ONE;

    sqP = a.f2Add (xP);                            // sqP = a + xP
    sqP = xP.f2Mult (sqP);                         // sqP = xP * sqP
    sqP = sqP.f2Add (F2elm.ONE);                   // sqP = sqP + 1
    sqQ = a.f2Add (xQ);                            // sqQ = a + xQ
    sqQ = xQ.f2Mult (sqQ);                         // sqQ = xq * sqQ
    sqQ = sqQ.f2Add (F2elm.ONE);                   // sqQ = sqQ + 1
    sqQ = xQ.f2Mult (sqQ);                         // sqQ = xq * sqQ
    sqP = xP.f2Mult (sqP);                         // sqP = xP * sqP
    sq = sqP.f2Mult (sqQ);                         // sq = sqP * sqQ

    // The first step after initial setup is ensure p and q have the right order

    for (j = 1; j < eB; j++) {
      numers = cubeIndeterminant (numers[0], numers[1], sq);
      denoms = cubeIndeterminant (denoms[0], denoms[1], sq);

      tlout = tplLine (p, q);
      p = new F2Point (tlout[0], tlout[1]);   q = new F2Point (tlout[2], tlout[3]);  

      alphan = uP.f2Mult (tlout[8]);             // alphan = alphan*uP
      alphan = uQD.f2Mult (alphan);              // alphan = alphan*uQD  
      alphad = uQ.f2Mult (tlout[9]);             // alphad = alphad*uQ 
      alphad = uPD.f2Mult (alphad);              // alphad = alphad*uPD
      t0 = uQD.f2Mult (uPD);                     // t0 = uQD*uPD
      betan = t0.f2Mult (tlout[10]);             // betan = betan*t0
      betad = t0.f2Mult (tlout[11]);             // betad = betad*t0
      uP = uP.f2Mult (tlout[4]);                 // uP = uP*xUP
      uPD = uPD.f2Mult (tlout[5]);               // uPD = uPD*zUP
      uQ = uQ.f2Mult (tlout[6]);                 // uQ = uQ*xUQ
      uQD = uQD.f2Mult (tlout[7]);               // uQD = uQD*zUQ

      numers = lineIndeterminantTpl (numers[0], numers[1], alphan, betan, sq);
      denoms = lineIndeterminantTpl (denoms[0], denoms[1], alphad, betad, sq);
    }

    numers = cubeIndeterminant (numers[0], numers[1], sq);
    denoms = cubeIndeterminant (denoms[0], denoms[1], sq);

    px = p.getX();  pz = p.getZ();
    qx = q.getX();  qz = q.getZ();

    t0 = a.f2Mult (pz);                     // t0 = a*pz
    t1 = px.f2Add (px);                     // t1 = px+px
    t1 = px.f2Add (t1);                     // t1 = t1+px
    t2 = t0.f2Add (t1);                     // t2 = t1+t0
    t1 = t2.f2Add (t0);                     // t1 = t2+t0
    lambdaP = t1.f2Mult(px);                // lambdaP = t1*px
    t1 = pz.f2Sqr();                        // t1 = pz^2
    lambdaP = lambdaP.f2Add (t1);           // lambdaP = lambdaP+t1
    t1 = t1.f2Mult (pz);                    // t1 = t1*pz
    t0 = sqP.f2Mult (t1);                   // t0 = t1*sqP
    t1 = t1.f2Mult (uPD);                   // t1 = t1*uPD
    t3 = uP.f2Sqr();                       // t3 = uP^2
    t0 = t0.f2Mult (t3);                    // t0 = t0*t3
    t0 = t0.f2Add (t0);                     // t0 = t0+t0
    t2 = t2.f2Mult (t0);                    // t2 = t2*t0
    t2 = t2.f2Add (t2);                     // t2 = t2+t2
    t3 = lambdaP.f2Mult (uPD);              // t3 = lambdaP*uPD
    t4 = t3.f2Sqr();                       // t4 = t3^2

    // Check order of p by asserting that 3^(eB-1)*p has order 3
    if ( !t2.f2Equals (t4) || t2.f2Equals (F2elm.ZERO))
      return false;

    t5 = a.f2Mult (qz);                     // t5 = a*qz 
    t6 = qx.f2Add (qx);                     // t6 = qx+qx  
    t6 = qx.f2Add (t6);                     // t6 = t6+qx
    t2 = t5.f2Add (t6);                     // t2 = t6+t5
    t6 = t2.f2Add (t5);                     // t6 = t2+t5
    lambdaQ = t6.f2Mult (qx);               // lambdaQ = t6*qx
    t6 = qz.f2Sqr();                       // t6 = qz^2
    lambdaQ = lambdaQ.f2Add (t6);           // lambdaQ = lambdaQ+t6
    t6 = t6.f2Mult (qz);                    // t6 = t6*qz
    t5 = sqQ.f2Mult (t6);                   // t5 = t6*sqQ
    t6 = t6.f2Mult (uQD);                   // t6 = t6*uQD
    t7 = uQ.f2Sqr();                       // t7 = uQ^2
    t5 = t5.f2Mult (t7);                    // t5 = t5*t7
    t5 = t5.f2Add (t5);                     // t5 = t5+t5
    t2 = t2.f2Mult (t5);                    // t2 = t2*t5
    t2 = t2.f2Add (t2);                     // t2 = t2+t2
    t7 = lambdaQ.f2Mult (uQD);              // t7 = lambdaQ*uQD
    t4 = t7.f2Sqr();                       // t4 = t7^2

    // Check order of q by asserting that 3^(eB-1)*q has order 3
    if ( !t2.f2Equals (t4) || t2.f2Equals (F2elm.ZERO))
      return false;

    // Next, check the Weil pairing to ensure that p and q are linearly independent

    lnQ = xQ.f2Mult (pz);                          // lnQ = xQ*pz
    lnQ = xP.f2Sub(lnQ);                           // lnQ = xP-lnQ
    lnQ = t3.f2Mult (lnQ);                         // lnQ = t3*lnQ
    lnQ = uPD.f2Mult(lnQ);                         // lnQ = lnQ*uPD
    lnQ = lnQ.f2Sub (t0);                          // lnQ = lnQ-t0
    ldP = xP.f2Mult (qz);                          // ldP = xP*qz
    ldP = xQ.f2Sub (ldP);                          // ldP = xQ-ldP
    ldP = t7.f2Mult (ldP);                         // ldP = t7*ldP
    ldP = uQD.f2Mult (ldP);                        // ldP = uQD*ldP
    ldP = ldP.f2Sub (t5);                          // ldP = ldP-t5
    lnP = uP.f2Mult (uQ);                          // lnP = uP*uQ
    lnP = lnP.f2Add (lnP);                         // lnP = lnP+lnP
    ldQ = sqP.f2Mult (lnP);                        // ldQ = lnP*sqP
    lnP = lnP.f2Mult (sqQ);                        // lnP = lnP*sqQ
    ldP = ldP.f2Mult (uP);                         // ldP = ldP*uP
    ldP = t1.f2Mult (ldP);                         // ldP = ldP*t1
    lnQ = lnQ.f2Mult (uQ);                         // lnQ = lnQ*uQ
    lnQ = t6.f2Mult (lnQ);                         // lnQ = lnQ*t6
    t1 = t1.f2Mult (t6);                           // t1 = t1*t6
    lnP = lnP.f2Mult (t1);                         // lnP = lnP*t1    
    ldQ = ldQ.f2Mult (t1);                         // ldQ = ldQ*t1
    t0 = numers[0];                                // t0 = alpha_numer
    numers[0] = lnP.f2Mult (t0);                   // alpha_numer = lnP*t0
    numers[0] = sqP.f2Mult (numers[0]);            // alpha_numer = alpha_numer*sqP
    t1 = lnQ.f2Mult (numers[1]);                   // t1 = lnQ*beta_numer
    numers[0] = numers[0].f2Add (t1);              // alpha_numer = t1+alpha_numer
    t1 = t0.f2Mult (sqQ);                          // t1 = t0*sqQ
    t1 = t1.f2Mult (lnQ);                          // t1 = t1*lnQ
    numers[1] = lnP.f2Mult (numers[1]);            // beta_numer = lnP*beta_numer
    numers[1] = t1.f2Add (numers[1]);              // beta_numer = beta_numer+t1
    t0 = denoms[0];                                // t0 = alpha_denom
    t1 = ldP.f2Mult (t0);                          // t1 = ldP*t0
    t1 = sqP.f2Mult (t1);                          // t1 = t1*sqP
    denoms[0] = denoms[1].f2Mult (ldQ);            // alpha_denom = ldQ*beta_denom
    denoms[0] = t1.f2Add (denoms[0]);              // alpha_denom = alpha_denom+t1
    t1 = t0.f2Mult (sqQ);                          // t1 = t0*sqQ
    t1 = ldQ.f2Mult (t1);                          // t1 = ldQ*t1
    denoms[1] = denoms[1].f2Mult (ldP);            // beta_denom = ldP*beta_denom
    denoms[1] = t1.f2Add (denoms[1]);              // beta_denom = beta_denom+t1    
    t2 = numers[0].f2Add (denoms[0]);              // t2 = alpha_numer+alpha_denom
    t2 = t2.f2Sqr();                              // t2 = t2^2
    t2 = t2.f2Mult (sqQ);                          // t2 = t2*sqQ
    t4 = numers[1].f2Add (denoms[1]);              // t4 = beta_numer+beta_denom
    t4 = t4.f2Sqr();                              // t4 = t4^2
    t4 = sqP.f2Mult (t4);                          // t4 = t4*sqP
    
    // Check if weil pairing = 1
    if (t2.f2Equals (t4))
      return false;

    // Next check that the third point is the difference between p and q

    t0 = xP.f2Add (xQ);                            // t0 = xP+xQ
    t1 = xQP.f2Mult (t0);                          // t1 = xQP*t0
    t1 = t1.f2Sub (F2elm.ONE);                     // t1 = t1-1
    t2 = xP.f2Mult (xQ);                           // t2 = xP*xQ
    t1 = t1.f2Add (t2);                            // t1 = t2+t1
    t1 = t1.f2Sqr();                              // t1 = t1^2
    t0 = t0.f2Add (xQP);                           // t0 = t0+xQP
    t0 = t0.f2Add (a);                             // t0 = t0+a
    t2 = t2.f2Mult (xQP);                          // t2 = t2*xQP
    t0 = t0.f2Mult (t2);                           // t0 = t0*t2
    t0 = t0.f2Add (t0);                            // t0 = t0+t0
    t0 = t0.f2Add (t0);                            // t0 = t0+t0

    if (!t0.f2Equals (t1))
      return false;

    return testCurve (eA, eB);
  }


  public boolean validateBsKey (int eA, int eB) {
    F2elm t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, fP, fQ, uP, uQ, vP, vQ;
    F2elm cP, cQ, alphaQi, betaPi, alphaPi, betaQi, alphaP, alphaQ, betaP, betaQ;
    F2elm a = pk.getA(), xP = pk.getP(), xQ = pk.getQ(), xQP = pk.getD(), px, pz, qx, qz;
    F2Point p, q;
    int j;
    MontCurve curve = new MontCurve (a, F2elm.ONE);

    px = xP;  pz = F2elm.ONE;
    qx = xQ;  qz = F2elm.ONE;

    fP = F2elm.ONE;
    fQ = F2elm.ONE;
    uP = F2elm.ONE;
    uQ = F2elm.ONE;
    vP = F2elm.ONE;
    vQ = F2elm.ONE;
    alphaP = F2elm.ZERO;
    alphaQ = F2elm.ZERO;
    betaP = F2elm.ONE;
    betaQ = F2elm.ONE;

    cQ = xQ.f2Add (a);                             // cQ = xQ+a
    cP = xP.f2Add (a);                             // cP = xP+a
    cQ = cQ.f2Mult(xQ);                            // cQ = cQ*xQ 
    cP = cP.f2Mult (xP);                           // cP = cP*xP 
    cQ = cQ.f2Add (F2elm.ONE);                     // cQ = cQ+1
    cP = cP.f2Add (F2elm.ONE);                     // cP = cP+1
    cQ = cQ.f2Mult (xQ);                           // cQ = cQ*xQ
    cP = cP.f2Mult (xP);                           // cP = cP*xP

    t0 = xP;
    t1 = xQ;

    for (j = 1; j < eA; j++) {
      t2 = px.f2Sqr();                          // t2 = px^2
      t11 = pz.f2Sqr();                         // t11 = pz^2
      t4 = qx.f2Sqr();                          // t4 = qx^2
      t10 = qz.f2Sqr();                         // t10 = qz^2
      t6 = t2.f2Sub (t11);                       // t6 = t2-t11
      betaPi = t2.f2Add (t2);                    // betaPi = t2+t2
      t2 = t2.f2Add (t11);                       // t2 = t2+t11
      t7 = t4.f2Sub (t10);                       // t7 = t4-t10
      alphaQi = t4.f2Add (t4);                   // alphaQi = t4+t4
      t4 = t4.f2Add (t10);                       // t4 = t4+t10
      t3 = px.f2Mult (pz);                       // t3 = px*pz 
      t5 = qx.f2Mult (qz);                       // t5 = qx*qz 
      t8 = a.f2Mult (t3);                        // t8 = a*t3
      t9 = a.f2Mult (t5);                        // t9 = a*t5
      t3 = t3.f2Add (t3);                        // t3 = t3+t3
      t5 = t5.f2Add (t5);                        // t5 = t5+t5
      betaPi = betaPi.f2Add (t8);                // betaPi = betaPi+t8
      t8 = t2.f2Add (t8);                        // t8 = t8+t2
      alphaQi = alphaQi.f2Add (t9);              // alphaQi = alphaQi+t9
      t9 = t4.f2Add (t9);                        // t9 = t9+t4
      t2 = a.f2Mult (t2);                        // t2 = a*t2
      t4 = a.f2Mult (t4);                        // t4 = a*t4
      betaPi = betaPi.f2Add (t8);                // betaPi = betaPi+t8
      alphaQi = alphaQi.f2Add (t9);              // alphaQi = alphaQi+t9
      betaPi = betaPi.f2Mult (t1);               // betaPi = betaPi*t1
      alphaQi = alphaQi.f2Mult (t0);             // alphaQi = alphaQi*t0
      t1 = px.f2Mult (t6);                       // t1 = px*t6
      t0 = qx.f2Mult (t7);                       // t0 = qx*t7
      t1 = t1.f2Sub (betaPi);                    // t1 = t1-betaPi
      t0 = t0.f2Sub (alphaQi);                   // t0 = t0-alphaQi
      betaPi = vP.f2Mult (t1);                   // betaPi = vP*t1
      alphaQi = vQ.f2Mult (t0);                  // alphaQi = vQ*t0
      t10 = qz.f2Mult (t10);                     // t10 = t10*qz
      t11 = pz.f2Mult (t11);                     // t11 = t11*pz
      t10 = t10.f2Mult (uQ);                     // t10 = t10*uQ
      t11 = t11.f2Mult (uP);                     // t11 = t11*uP
      betaPi = betaPi.f2Mult (t10);              // betaPi = betaPi*t10
      alphaQi = alphaQi.f2Mult (t11);            // alphaQi = alphaQi*t11
      t0 = t10.f2Mult (t11);                     // t10 = t10*t11
      t10 = t10.f2Add (t10);                     // t10 = t10+t10
      alphaPi = cQ.f2Mult (t10);                 // alphaPi = cQ*t10
      betaQi = cP.f2Mult (t10);                  // betaQi = cP*t10
      uQ = uQ.f2Mult (t7);                       // uQ = uQ*t7
      uP = uP.f2Mult (t6);                       // uP = uP*t6
      t8 = t8.f2Add (t8);                        // t8 = t8+t8
      t9 = t9.f2Add (t9);                        // t9 = t9+t9
      pz = t3.f2Mult (t8);                       // pz = t3*t8
      qz = t5.f2Mult (t9);                       // qz = t5*t9
      t8 = t8.f2Mult (px);                       // t8 = t8*px
      t9 = t9.f2Mult (qx);                       // t9 = t9*qx
      px = t6.f2Sqr();                          // px = t6^2
      qx = t7.f2Sqr();                          // qx = t7^2
      t4 = t4.f2Add (t5);                        // t4 = t4+t5
      t4 = t5.f2Add (t4);                        // t4 = t4+t5
      t2 = t2.f2Add (t3);                        // t2 = t2+t3
      t2 = t3.f2Add (t2);                        // t2 = t2+t3
      t4 = t4.f2Mult (t5);                       // t4 = t4*t5
      t2 = t2.f2Mult (t3);                       // t2 = t2*t3
      t4 = t4.f2Add (qx);                        // t4 = t4+qx
      t2 = t2.f2Add (px);                        // t2 = t2+px
      uQ = uQ.f2Mult (t4);                       // uQ = uQ*t4
      uP = uP.f2Mult (t2);                       // uP = uP*t2
      t9 = t9.f2Sqr();                          // t9 = t9^2
      t8 = t8.f2Sqr();                          // t8 = t8^2
      vQ = vQ.f2Mult (t9);                       // vQ = vQ*t9
      vP = vP.f2Mult (t8);                       // vP = vP*t8
      vQ = vQ.f2Add (vQ);                        // vQ = vQ+vQ
      vP = vP.f2Add (vP);                        // vP = vP+vP
      t4 = alphaP.f2Sqr();                      // t4 = alphaP^2
      t5 = betaP.f2Sqr();                       // t5 = betaP^2
      t6 = alphaP.f2Mult (betaP);                // t6 = alphaP*betaP
      t6 = t6.f2Add (t6);                        // t6 = t6+t6 
      t4 = t4.f2Mult (cP);                       // t4 = t4*cP
      t5 = t5.f2Mult (cQ);                       // t5 = t5*cQ
      t4 = t4.f2Add (t5);                        // t4 = t4+t5
      alphaP = alphaPi.f2Mult (t4);              // alphaP = alphaPi*t4
      betaP = betaPi.f2Mult (t4);                // betaP = betaPi*t4
      t4 = betaPi.f2Mult (t6);                   // t4 = t6*betaPi
      t4 = t4.f2Mult (cQ);                       // t4 = t4*cQ
      t6 = t6.f2Mult (alphaPi);                  // t6 = t6*alphaPi
      t6 = cP.f2Mult (t6);                       // t6 = t6*cP
      alphaP = alphaP.f2Add (t4);                // alphaP = alphaP+t4
      betaP = betaP.f2Add (t6);                  // betaP = betaP+t6
      t4 = alphaQ.f2Sqr();                      // t4 = alphaQ^2
      t5 = betaQ.f2Sqr();                       // t5 = betaQ^2
      t6 = alphaQ.f2Mult (betaQ);                // t6 = alphaQ*betaQ
      t6 = t6.f2Add (t6);                        // t6 = t6+t6
      t4 = t4.f2Mult (cP);                       // t4 = t4*cP
      t5 = t5.f2Mult (cQ);                       // t5 = t5*cQ
      t4 = t4.f2Add (t5);                        // t4 = t4+t5
      alphaQ = alphaQi.f2Mult (t4);              // alphaQ = alphaQi*t4
      betaQ = betaQi.f2Mult (t4);                // betaQ = betaQi*t4
      t4 = betaPi.f2Mult (t6);                   // t4 = t6*betaPi
      t4 = cQ.f2Mult (t4);                       // t4 = t4*cQ
      t5 = t6.f2Mult (betaQi);                   // t5 = t6*betaQi
      t5 = t5.f2Mult (cQ);                       // t5 = t5*cQ
      alphaQ = alphaQ.f2Add (t5);                // alphaQ = alphaQ+t5
      t5 = t6.f2Mult (alphaQi);                  // t5 = t6*alphaQi
      t5 = cP.f2Mult (t5);                       // t5 = t5*cP
      betaQ = betaQ.f2Add (t5);                  // betaQ = betaQ+t5
      t0 = xP.f2Mult (qz);                       // t0 = xP*qz 
      t1 = xQ.f2Mult (pz);                       // t1 = xQ*pz 
      t2 = qx.f2Sub (t0);                        // t2 = qx-t0
      t3 = px.f2Sub (t1);                        // t3 = px-t1
      t2 = t2.f2Mult (pz);                       // t2 = t2*pz
      t3 = t3.f2Mult (qz);                       // t3 = t3*qz
      alphaP = alphaP.f2Mult (t2);               // alphaP = alphaP*t2
      betaP = betaP.f2Mult (t2);                 // betaP = betaP*t2
      alphaQ = alphaQ.f2Mult (t3);               // alphaQ = alphaQ*t3
      betaQ = betaQ.f2Mult (t3);                 // betaQ = betaQ*t3
    }

    t2 = xQ.f2Mult (pz);                           // t2 = xQ*pz
    t3 = xP.f2Mult (qz);                           // t3 = xP*ZQ
    t2 = px.f2Sub (t2);                            // t2 = px-t2
    t3 = qx.f2Sub (t3);                            // t3 = qx-t3
    t2 = t2.f2Mult (qz);                           // t2 = t2*qz
    t3 = t3.f2Mult (pz);                           // t3 = t3*pz
    t4 = alphaP.f2Sqr();                          // t4 = alphaP^2
    t5 = betaP.f2Sqr();                           // t5 = betaP^2
    t6 = alphaP.f2Mult (betaP);                    // t6 = alphaP*betaP
    t6 = t6.f2Add (t6);                            // t6 = t6+t6
    t7 = alphaQ.f2Sqr();                          // t7 = alphaQ^2
    t8 = betaQ.f2Sqr();                           // t8 = betaQ^2
    t9 = alphaQ.f2Mult (betaQ);                    // t9 = alphaQ*betaQ
    t9 = t9.f2Add (t9);                            // t9 = t9+t9
    t4 = t4.f2Mult (cP);                           // t4 = t4*cP
    t5 = t5.f2Mult (cQ);                           // t5 = t5*cQ
    t7 = t7.f2Mult (cP);                           // t7 = t7*cP
    t8 = t8.f2Mult (cQ);                           // t8 = t8*cQ
    t4 = t4.f2Add (t5);                            // t4 = t4+t5
    t7 = t7.f2Add (t8);                            // t7 = t7+t8
    t4 = t2.f2Mult (t4);                           // t4 = t2*t4
    t7 = t3.f2Mult (t7);                           // t7 = t3*t7
    t7 = t4.f2Sub (t7);                            // t7 = t4-t7
    t7 = t7.f2Sqr();                              // t7 = t7^2
    t3 = t3.f2Mult (t9);                           // t3 = t3*t9
    t2 = t2.f2Mult (t6);                           // t2 = t2*t6
    t3 = t3.f2Sub (t2);                            // t3 = t3-t2
    t3 = t3.f2Sqr();                              // t3 = t3^2
    t3 = cP.f2Mult (t3);                           // t3 = t3*cP
    t3 = cQ.f2Mult (t3);                           // t3 = t3*cQ
    
    // Check the order of q by asserting that q^(eA-1) has order 2
    q = new F2Point (qx, qz);
    q = curve.xDbl (q);
    if (qz.f2Equals (F2elm.ZERO) || !(q.getZ().f2Equals (F2elm.ZERO)))
        return false;

    // Check the order of p by asserting that q^(eA-1) has order 2
    p = new F2Point (px, pz);
    p = curve.xDbl (p);
    if (qz.f2Equals (F2elm.ZERO) || !(q.getZ().f2Equals (F2elm.ZERO)))
      return false;

    // Check that the Weil pairing is nontrivial to ensure p and q are independent
    if (t3.f2Equals (t7))
      return false;

    // Check that the third point is the difference of p and q
    t0 = xP.f2Add (xQ);                            // t0 = xP+xQ
    t1 = xQP.f2Mult (t0);                          // t1 = xQP*t0
    t1 = t1.f2Sub (F2elm.ONE);                     // t1 = t1-1
    t2 = xP.f2Mult (xQ);                           // t2 = xP*xQ
    t1 = t1.f2Add (t2);                            // t1 = t2+t1
    t1 = t1.f2Sqr();                              // t1 = t1^2
    t0 = t0.f2Add (xQP);                           // t0 = t0+xQP
    t0 = a.f2Add (t0);                             // t0 = t0+a
    t2 = xQP.f2Mult (t2);                          // t2 = t2*xQP
    t0 = t0.f2Mult (t2);                           // t0 = t0*t2
    t0 = t0.f2Add (t0);                            // t0 = t0+t0
    t0 = t0.f2Add (t0);                            // t0 = t0+t0

    if (!t0.f2Equals (t1))
      return false;

    return testCurve (eA, eB);
  }
}



