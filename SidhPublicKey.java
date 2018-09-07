package sidh;

/**************************************************************************************************
 *
 * Implements public keys for the Supersingular Isogeny Diffie-Hellman key exchange algorithm. 
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;

class SidhPublicKey {
  F2elm phiPx;
  F2elm phiQx;
  F2elm phiDx;


  public SidhPublicKey (int aOrB, SidhPrivateKey k, SidhKeyExchange params) {
    if (aOrB == SidhKeyExchange.PARTYA) 
      genPubKeyA (k, params);
    else 
      genPubKeyB (k, params);
  } 


  public SidhPublicKey (SidhPublicKey k) {
    phiPx = new F2elm (k.phiPx);
    phiQx = new F2elm (k.phiQx);
    phiDx = new F2elm (k.phiDx);
  }


  public SidhPublicKey (byte[] inBytes) {
    int len = (inBytes.length) / 3;
    phiPx = new F2elm (Arrays.copyOfRange (inBytes, 0, len));
    phiQx = new F2elm (Arrays.copyOfRange (inBytes, len, 2*len));
    phiDx = new F2elm (Arrays.copyOfRange (inBytes, 2*len, 3*len));
  }


  public F2elm getP () {
    return phiPx;
  }


  public F2elm getQ () {
    return phiQx;
  }


  public F2elm getD () {
    return phiDx;
  }
    
    
  public byte[] serialize() {
    int f2size = 2 * Felm.primesize;
    byte[] retval = new byte[3*f2size];
    
    System.arraycopy (phiPx.toByteArray(), 0, retval, 0, f2size);
    System.arraycopy (phiQx.toByteArray(), 0, retval, f2size, f2size);
    System.arraycopy (phiDx.toByteArray(), 0, retval, 2*f2size, f2size);

    return retval;
  }


  public boolean publicKeyEquals (SidhPublicKey k2) {
    if (phiPx.f2Equals (k2.phiPx) == false)
      return false;
    if (phiQx.f2Equals (k2.phiQx) == false)
      return false;
    return phiDx.f2Equals (k2.phiDx);
  }
    

  public int hashCode() {
    return Arrays.hashCode (serialize());
  }


  protected void genPubKeyA (SidhPrivateKey privKey, SidhKeyExchange params) {
    // Given A's private key compute the corresponding public key

    MontCurve curve;
    FourIsogeny fourIsog;
    F2Point r, phiP, phiQ, phiD, pts[];
    F2elm invs[], coeffs[], genA[], genB[];
    int maxIntPointsA, maxA, splitsA[], obits, row, index = 0, npts = 0, ptsIdx[], m, i, ii = 0;
 
    obits = params.getObitsA();
    curve = new MontCurve (params.getCurve());
    curve.updateA24();
    
    genA = params.getGenA();
    genB = params.getGenB();

    r = curve.ladder3pt (genA[0], genA[1], genA[2], privKey.getKey(), obits);
    
    phiP = new F2Point (genB[0], F2elm.ONE);
    phiQ = new F2Point (genB[1], F2elm.ONE);
    phiD = new F2Point (genB[2], F2elm.ONE);
    
    fourIsog = new FourIsogeny (curve);
    fourIsog.setA24plus (new F2elm (F2elm.ONE));
    fourIsog.setC24 (F2elm.leftShift (F2elm.ONE, 1));

    maxIntPointsA = params.getMIPA();
    maxA = params.getMaxA();
    splitsA = params.getSplitsA();
    pts = new F2Point[maxIntPointsA];
    ptsIdx = new int[maxIntPointsA];
    
    for (row = 1; row < maxA; row++) {
      while (index < maxA - row) {
	pts[npts] = r;
        ptsIdx[npts++] = index;
        m = splitsA[ii++];
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

    invs = F2elm.inv3Way (phiP.getZ(), phiQ.getZ(), phiD.getZ());

    phiPx = F2elm.mult (invs[0], phiP.getX());
    phiQx = F2elm.mult (invs[1], phiQ.getX());
    phiDx = F2elm.mult (invs[2], phiD.getX());
  }


  protected void genPubKeyB (SidhPrivateKey privKey, SidhKeyExchange params) {
    // Given B's private key, compute the corresponding public key

    MontCurve curve; 
    F2elm genA[], genB[], invs[];
    int maxIntPointsB, maxB, splitsB[], obits, row, m, index = 0, ptsIdx[], npts = 0, i, ii = 0;
    
    F2Point r, phiP, phiQ, phiD, pts[];
    ThreeIsogeny threeIsog;

    obits = params.getObitsB();
    curve = new MontCurve (params.getCurve());
    curve.updateA24();

    genA = params.getGenA();
    genB = params.getGenB();
    
    r = curve.ladder3pt (genB[0], genB[1], genB[2], privKey.getKey(), obits);

    phiP = new F2Point (genA[0], F2elm.ONE);
    phiQ = new F2Point (genA[1], F2elm.ONE);
    phiD = new F2Point (genA[2], F2elm.ONE);
    
    threeIsog = new ThreeIsogeny (curve);
    threeIsog.setA24plus (F2elm.leftShift (F2elm.ONE, 1));
    threeIsog.setA24minus (F2elm.negate (threeIsog.getA24plus()));

    maxIntPointsB = params.getMIPB();
    maxB = params.getMaxB();
    splitsB = params.getSplitsB();
    pts = new F2Point[maxIntPointsB];
    ptsIdx = new int[maxIntPointsB];
    
    for (row = 1; row < maxB; row++) {
      while (index < maxB - row) {
        pts[npts] = r;
        ptsIdx[npts++] = index;
        m = splitsB[ii++];
        r = threeIsog.xTple (r, m);
        index += m;
      }

      threeIsog.get3Isog (r);

      for (i = 0; i < npts; i++) 
        pts[i] = threeIsog.eval3Isog (pts[i]);
      
      phiP = threeIsog.eval3Isog (phiP);      
      phiQ = threeIsog.eval3Isog (phiQ);
      phiD = threeIsog.eval3Isog (phiD);

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    threeIsog.get3Isog (r);

    phiP = threeIsog.eval3Isog (phiP);
    phiQ = threeIsog.eval3Isog (phiQ);
    phiD = threeIsog.eval3Isog (phiD);
    
    invs = F2elm.inv3Way (phiP.getZ(), phiQ.getZ(), phiD.getZ());

    phiPx = F2elm.mult (invs[0], phiP.getX());
    phiQx = F2elm.mult (invs[1], phiQ.getX());
    phiDx = F2elm.mult (invs[2], phiD.getX());
  }
}
