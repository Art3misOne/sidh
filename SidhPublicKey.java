package sidh;

/**************************************************************************************************
 *
 * Implements public keys for the Supersingular Isogeny Diffie-Hellman key exchange algorithm. 
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;

class SidhPublicKey {
  F2elm a;
  F2elm phiPx;
  F2elm phiQx;
  F2elm phiDx;


  public SidhPublicKey (int aOrB, SidhPrivateKey k, SidhKeyExchange params) {
    if (aOrB == SidhKeyExchange.PARTYA)
      genPubKeyA (k, params);
    else
      genPubKeyB (k, params);
  } 


  public SidhPublicKey (F2elm ain, F2elm phip, F2elm phiq, F2elm phiqp) {
    a = ain;
    phiPx = phip;
    phiQx = phiq;
    phiDx = phiqp;
  }


  public SidhPublicKey (SidhPublicKey k) {
    a = k.a;
    phiPx = k.phiPx;
    phiQx = k.phiQx;
    phiDx = k.phiDx;
  }


  public SidhPublicKey (byte[] inBytes) {
    int len = (inBytes.length) / 4;
    a = new F2elm (Arrays.copyOfRange (inBytes, 0, len));
    phiPx = new F2elm (Arrays.copyOfRange (inBytes, len, 2*len));
    phiQx = new F2elm (Arrays.copyOfRange (inBytes, 2*len, 3*len));
    phiDx = new F2elm (Arrays.copyOfRange (inBytes, 3*len, 4*len));
  }


  public F2elm getA() {
    return a;
  }


  public F2elm getP() {
    return phiPx;
  }


  public F2elm getQ() {
    return phiQx;
  }


  public F2elm getD() {
    return phiDx;
  }


  public byte[] serialize() {
    int f2size = 2 * Felm.primesize;
    byte[] retval = new byte[4*f2size];
    System.arraycopy (a.toByteArray(), 0, retval, 0, f2size);
    System.arraycopy (phiPx.toByteArray(), 0, retval, f2size, f2size);
    System.arraycopy (phiQx.toByteArray(), 0, retval, 2*f2size, f2size);
    System.arraycopy (phiDx.toByteArray(), 0, retval, 3*f2size, f2size);
    return retval;
  }


  public boolean equals (SidhPublicKey pk) {  
    if (a.f2Equals (pk.a))
      if (phiPx.f2Equals (pk.phiPx))
        if (phiQx.f2Equals (pk.phiQx))
	    return (phiDx.f2Equals (pk.phiDx));
    return false;
  }


  public int hashCode() {
    return Arrays.hashCode (serialize());
  }


  private static F2Point secretPt (FpPointAffine p, SidhPrivateKey pk, MontCurve curve, int obits) {
    // Computes key generation entirely in the base field by exploiting a 1-dimensional 
    //   Montgomery ladder in the trace 0 subgroup and recovering the y-coordinate for the 
    //   addition.
    // Inputs: Scalar m (secret key for A or B), affine point p = (x,y), (let q = (-x,y)) 
    // Outputs: r = ( (rx0+rx1*i) : rz0 ) where rx/rz = the x coordinate of p + mq

    int nbits;
    FpPointProjective st[];
    Felm rx0, rx1, rz0, rz1, sx, sz, tx, tz, x, y, x1, y1, t0, t1, t2;

    x = p.getX();
    y = p.getY();
    x1 = x.fpNegate();
    y1 = y;

    st = new FpPointProjective[2];

    st = curve.ladder (x1, pk.getKey(), obits);

    sx = st[0].getX();
    sz = st[0].getZ();
    tx = st[1].getX();
    tz = st[1].getZ();

    // rx0 := 2*y*y1*sz^2*tz - tz*(sx*x1+sz)*(sx+x1*sz) + tx*(sx-x1*sz)^2 
    //       - 4*y1^2*sz*tz^2*(sx+x*sz)*(sx-x*sz)^2;
    // rx1 := 4*y*y1*sz^2*tz*(tz*(sx*x1+sz)*(sx+x1*sz) - tx*(sx-x1*sz)^2);
    // rz0 := 4*y1^2*sz^2*tz^2*(sx-x*sz)^2;

    rx1 = x1.fpMult(sz);                     // rx1 = x1 * sz
    rx0 = sx.fpMult(x1);                     // rx0 = sx * x1
    t0 = sx.fpSub(rx1);                      // t0 = sx - (x1*sz)
    rx1 = sx.fpAdd(rx1);                     // rx1 = sx + (x1*sz)
    t0 = t0.fpSqr();                         // t0 = (sx-x1*sz) ^ 2
    rx0 = rx0.fpAdd(sz);                     // rx0 = sx*x1 + sz
    t0 = t0.fpMult(tx);                      // t0 = tx * (sx-x1*sz)^2
    rx0 = rx0.fpMult(rx1);                   // rx0 = (sx*x1+sz) * (sx+x1*sz) 
    t2 = y1.fpMult(tz);                      // t2 = y1 * tz
    t1 = y.fpMult(sz);                       // t1 = y * sz
    t2 = t2.fpAdd(t2);                       // t2 = 2 * y1*tz
    rx1 = t2.fpMult(sz);                     // rx1 = 2*y1*tz * sz
    rx0 = rx0.fpMult(tz);                    // rx0 = (sx*x1+sz)*(sx+x1*sz) * tz 
    rx0 = rx0.fpSub(t0);                     // rx0 = (sx*x1+sz)*(sx+x1*sz)*tz - tx*(sx-x1*sz)^2  
    t1 = t1.fpMult(rx1);                     // t1 = y*sz *  2*y1*tz*sz = 2*y*y1*sz^2*tz
    t0 = rx1.fpSqr();                        // t0 = (2*y1*tz*sz) ^ 2
    t2 = t2.fpMult(rx1);                     // t2 = 2*y1*tz * 2*y1*tz*sz = 4*y1^2*sz*tz^2
    rx1 = t1.fpMult(rx0);                    // rx1 = 2*y*y1*sz^2*tz * 
                                             //       (tz*(sx*x1+sz)*(sx+x1*sz) - tx*(sx-x1*sz)^2)
    rz0 = t1.fpAdd(rx0);                     // rz0 = 2*y*y1*sz^2*tz + (sx*x1+sz)*(sx+x1*sz)*tz - 
                                             //       tx*(sx-x1*sz)^2  
    rx1 = rx1.fpAdd(rx1);                    // rx1 = 2 * 2*y*y1*sz^2*tz*(tz*(sx*x1+sz)*(sx+x1*sz) 
                                             //       - tx*(sx-x1*sz)^2) 
                                             //     matching rx1 above
    t1 = t1.fpSub(rx0);                      // t1 = 2*y*y1*sz^2*tz - tz*(sx*x1+sz)*(sx+x1*sz) + 
                                             //      tx*(sx-x1*sz)^2
    rx0 = x.fpMult(sz);                      // rx0 = x * sz
    t1 = t1.fpMult(rz0);                     // t1 = (2*y*y1*sz^2*tz - (sx*x1+sz)*(sx+x1*sz)*tz + 
                                             //      tx*(sx-x1*sz)^2)  
    rz0 = sx.fpSub(rx0);                     // rz0 = sx - x*sz
    rx0 = sx.fpAdd(rx0);                     // rx0 = sx + x*sz
    rz0 = rz0.fpSqr();                       // rz0 = (sx-x*sz) ^ 2
    t2 = t2.fpMult(rx0);                     // t2 = 4*y1^2*sz*tz^2 * (sx+x*sz) 
    t2 = t2.fpMult(rz0);                     // t2 = 4*y1^2*sz*tz^2*(sx+x*sz) * (sx-x*sz)^2  
    rz0 = rz0.fpMult(t0);                    // rz0 = (sx-x*sz)^2 * (2*y1*tz*sz)^2
                                             //     = 4*y1^2*sz^2*tz^2*(sx-x*sz)^2 
                                             //     matching rz0 above
    rx0 = t1.fpSub(t2);                      // rx0 = 2*y*y1*sz^2*tz - tz*(sx*x1+sz)*(sx+x1*sz) + 
                                             //       tx*(sx-x1*sz)^2) - 
                                             //       4*y1^2*sz*tz^2*(sx+x*sz)*(sx-x*sz)^2  
                                             //     matching rx0 above
    rz1 = Felm.ZERO;

    return new F2Point (rx0, rx1, rz0, rz1);
  }


  protected void genPubKeyA (SidhPrivateKey privKeyA, SidhKeyExchange params) {
    // Given A's private key and the relevant key exchange parameters, compute her 
    // corresponding public key.

    MontCurve curve; 
    FpPointAffine generatorA, generatorB; 
    int maxIntPointsA, maxA, splitsA[], obits;
 
    // Extract key exchange parameters
    curve = params.getCurve();
    generatorA = params.getGeneratorA();
    generatorB = params.getGeneratorB();
    maxIntPointsA = params.getMIPA();
    maxA = params.getMaxA();
    splitsA = params.getSplitsA();
    obits = params.getObits();

    int row, index = 0, npts = 0, ptsIdx[], m, i;
    F2Point r, phiP, phiQ, phiD, pts[];
    F2elm invs[], coeffs[];
    FourIsogeny fourIsog = new FourIsogeny (curve);

    r = secretPt (generatorA, privKeyA, curve, obits);

    // Set phiP = (pBx : 1)
    phiPx = new F2elm (generatorB.getX(), Felm.ZERO);
    phiP = new F2Point (phiPx, F2elm.ONE);

    // Set phiQ = (-pBx : 1)
    phiQx = phiPx.f2Negate();
    phiQ = new F2Point (phiQx, F2elm.ONE);

    // Set phiD = QB - PB
    phiD = MontCurve.fpDistortAndDiff (phiPx.f2Get0());

    phiP = fourIsog.first4Isog (phiP);
    fourIsog = new FourIsogeny (curve);      // reset curve parameters
    phiQ = fourIsog.first4Isog (phiQ);
    fourIsog = new FourIsogeny (curve);      // reset curve parameters
    phiD = fourIsog.first4Isog (phiD);
    fourIsog = new FourIsogeny (curve);      // reset curve parameters
    r = fourIsog.first4Isog (r);

    pts = new F2Point[maxIntPointsA];
    ptsIdx = new int[maxIntPointsA];

    for (row = 1; row < maxA; row++) {
      while (index < maxA - row) {
        pts[npts] = r;
        ptsIdx[npts] = index;
        npts++;
        m = splitsA[maxA - index - row];
        r = fourIsog.xDble (r, 2*m);
        index += m;
      }

      fourIsog.get4Isog (r);

      for (i = 0; i < npts; i++)
        pts[i] = fourIsog.eval4Isog (pts[i]);

      phiP = fourIsog.eval4Isog (phiP);
      phiQ = fourIsog.eval4Isog (phiQ);
      phiD = fourIsog.eval4Isog (phiD);

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }
       
    fourIsog.get4Isog (r);

    phiP = fourIsog.eval4Isog (phiP);
    phiQ = fourIsog.eval4Isog (phiQ);
    phiD = fourIsog.eval4Isog (phiD);

    invs = F2elm.inv4Way (fourIsog.getC(), phiP.getZ(), phiQ.getZ(), phiD.getZ());

    a = invs[0].f2Mult (fourIsog.getA());    
    phiPx = invs[1].f2Mult (phiP.getX());
    phiQx = invs[2].f2Mult (phiQ.getX());
    phiDx = invs[3].f2Mult (phiD.getX());
  }


  protected void genPubKeyB (SidhPrivateKey privKeyB, SidhKeyExchange params) {
    // Given B's private key, compute the corresponding public key

    MontCurve curve; 
    FpPointAffine generatorA, generatorB; 
    int maxIntPointsB, maxB, splitsB[], obits;
 
    // Extract key exchange parameters
    curve = params.getCurve();
    generatorA = params.getGeneratorA();
    generatorB = params.getGeneratorB();
    maxIntPointsB = params.getMIPB();
    maxB = params.getMaxB();
    splitsB = params.getSplitsB();
    obits = params.getObits();

    F2Point r, phiP, phiQ, phiD, pts[];
    int i, row, m, index = 0, ptsIdx[], npts = 0;
    F2elm invs[];
    ThreeIsogeny threeIsog = new ThreeIsogeny (curve);

    pts = new F2Point[maxIntPointsB];
    ptsIdx = new int[maxIntPointsB];

    r = secretPt (generatorB, privKeyB, curve, obits);

    // Set phiP = (pAx : 1)
    phiPx = new F2elm (generatorA.getX(), Felm.ZERO);
    phiP = new F2Point (phiPx, F2elm.ONE);

    // Set phiQ = (-pAx : 1)
    phiQx = phiPx.f2Negate();
    phiQ = new F2Point (phiQx, F2elm.ONE);

    // Set phiD = QB - PB
    phiD = MontCurve.fpDistortAndDiff (phiPx.f2Get0());

    for (row = 1; row < maxB; row++) {
      while (index < maxB - row) {
        pts[npts] = r;
        ptsIdx[npts] = index;
        npts++;
        m = splitsB[maxB - index - row];
        r = threeIsog.xTple (r, m);
        index += m;
      }

      threeIsog.get3Isog (r);

      for (i = 0; i < npts; i++) 
        pts[i] = threeIsog.eval3Isog (r, pts[i]);

      phiP = threeIsog.eval3Isog (r, phiP);
      phiQ = threeIsog.eval3Isog (r, phiQ);
      phiD = threeIsog.eval3Isog (r, phiD);

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    threeIsog.get3Isog (r);

    phiP = threeIsog.eval3Isog (r, phiP);
    phiQ = threeIsog.eval3Isog (r, phiQ);
    phiD = threeIsog.eval3Isog (r, phiD);

    invs = F2elm.inv4Way (threeIsog.getC(), phiP.getZ(), phiQ.getZ(), phiD.getZ());

    a = invs[0].f2Mult (threeIsog.getA());    
    phiPx = invs[1].f2Mult (phiP.getX());
    phiQx = invs[2].f2Mult (phiQ.getX());
    phiDx = invs[3].f2Mult (phiD.getX());
  }
}
