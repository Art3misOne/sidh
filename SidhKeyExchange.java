
package sidh;

/**************************************************************************************************
 *
 * Implements supersingular isogeny Diffie-Hellman key exchange. Two constructors are provided. If
 * no arguments are passed, then it will default to using the parameters recommended by Microsoft
 * Research in the paper "Efficient algorithms for supersingular isogeny Diffie-Hellman" by 
 * Costello, Longa, and Naerhig. Alternate key exchange parameters can be used by passing them as
 * arguments to the constructor.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.lang.System;

class SidhKeyExchange {
  public static final int PARTYA = 0;
  public static final int PARTYB = 1;

  String isogenyId;

  int f;
  int lA;
  int lB;
  int eA;
  int eB;

  BigInteger prime;
  BigInteger orderA;
  BigInteger orderB;

  int obits;

  FpPointAffine pA;
  FpPointAffine pB;

  int maxA;
  int maxB;

  int maxIntPointsA;
  int maxIntPointsB;

  int[] splitsA;
  int[] splitsB;

  MontCurve baseCurve;

  public SidhKeyExchange (String id, int fIn, int laIn, int lbIn, int eaIn, int ebIn, 
			  FpPointAffine paIn, FpPointAffine pbIn, int mA, int mB, int mIPA, 
                          int mIPB, int[] sA, int[] sB, MontCurve curve) {                       
    int i;
    BigInteger temp;

    isogenyId = id;

    f = fIn;
    lA = laIn;
    lB = lbIn;
    eA = eaIn;
    eB = ebIn;

    orderA = BigInteger.valueOf(lA);     
    orderA = (orderA).pow(eA);               // orderA = lA^eA

    orderB = BigInteger.valueOf(lB);
    orderB = orderB.pow(eB);                 // orderB = lB^eB

    // obits = smallest multiple of 32 larger than number of bits in max (orderA, orderB)
    temp = orderA.max(orderB);
    obits = (temp.bitLength()-1) / 32;
    obits = (obits * 32) + 1;

    temp  = BigInteger.valueOf(f);
    prime = orderA.multiply(orderB);
    prime = prime.multiply(temp);
    prime = prime.subtract(BigInteger.ONE);  // prime = f * lA^eA * lB^eB - 1

    try {
      Felm.setPrime (prime);
    } catch (InvalidFieldException ex) {
      System.out.println ("\nParameters result in an invalid field. Using default parameters.");
      setDefaultParameters();
      return;
    }

    pA = paIn;
    pB = pbIn;

    maxA = mA;
    maxB = mB;
    
    maxIntPointsA = mIPA;
    maxIntPointsB = mIPB;

    splitsA = new int[mA];
    splitsB = new int[mB];

    System.arraycopy (sA, 0, splitsA, 0, mA);
    System.arraycopy (sB, 0, splitsB, 0, mB);

    baseCurve = new MontCurve (curve);
  }


  public SidhKeyExchange() {
    setDefaultParameters();
  }


  public void setDefaultParameters() {
    // Default parameters from "Efficient Algorithms for Supersingular Isogeny Diffie-Hellman" 
    // by Costello, Longa, and Naehrig

    BigInteger pAx, pAy, pBx, pBy;

    isogenyId = "sidhp751";

    f = 1;
    lA = 2;
    lB = 3;
    eA = 372;
    eB = 239;

    obits = 384;                             // Smallest multiple of 32 greater than order length

    // prime = f*(lA^eA)*(lB^eB) - 1
    prime = new BigInteger("6fe5d541f71c0e12909f97badc668562b5045cb25748084e9867d6ebe876da9" +
                           "59b1a13f7cc76e3ec968549f878a8eeafffffffffffffffffffffffffffffffff" +
                           "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);

    try {
      Felm.setPrime (prime);
    } catch (InvalidFieldException ex) {
      System.out.println ("\nUnrecoverable error: Change in default parameters resulted in an " +
                            "invalid field.\n");
      System.exit (-99); 
    }

    // orderA = lA^eA, orderB = lB^eB
    orderA = new BigInteger ("1000000000000000000000000000000000000000000000000000000000000" +
                             "000000000000000000000000000000000", 16);
    orderB = new BigInteger ("6fe5d541f71c0e12909f97badc668562b5045cb25748084e9867d6ebe876d" +
                             "a959b1a13f7cc76e3ec968549f878a8eeb", 16);

    pAx = new BigInteger ("3e82027a38e9429c8d36ff46bcc93fa23f89f6be06d2b1317ad90438621783fd" +
                          "b7a4ad3e83e86cae096d5db822c98e561e008fa0e3f3b9ac2f40c56d6fa4a58a20" +
                          "449af1f1335661d14ab7347693632646086ce3acd54b0346f5cce233e9", 16);
    pAy = new BigInteger ("3bbf8dcd4e7eb6236f5f598d56eb5e15915a755883b7c331b043da010e6a163a" +
                          "7421dfa8378d1e911f50bf3f721a8ed5950d80325a8d0f147ef3bd0cfec5236c7f" +
                          "ac9e69f7fd5a99ebec3b5b8b000f8eea737089343012e0d620bfb341d5", 16);
    pBx = new BigInteger ("2f1d80ef06ef960a01ab8ff409a2f8d5bce859ed725de145fe2d525160e0a3ad" +
                          "8e17b9f9238cd5e69cf26df237429bd3778659023b9ecb610e30288a7770d3785a" +
                          "aaa4d646c576aecb94b919aeedd9e1df566c1d26d376ed2325dcc93103", 16);
    pBy = new BigInteger ("127a46d082a1acaf351f09ab55a15445287ed1cc55dc35892123951d4b6e302c" +
                          "5129c049eeb399a6edb2eeb2f9b0a94f06cdfb3eade76eba0c8419745e97d12754" +
                          "f00e898a315b529122cfe3ca6bbc6baf5f6ba40bb91479226a0687894", 16);

    pA = new FpPointAffine (pAx, pAy);
    pB = new FpPointAffine (pBx, pBy);

    maxA = 185;
    maxB = 239;

    maxIntPointsA = 8;
    maxIntPointsB = 10;

    splitsA = new int[] { 
        0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 10, 9, 12, 
        11, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 19, 19, 
        17, 18, 19, 20, 21, 22, 21, 23, 22, 24, 24, 25, 25, 27, 27, 27, 28, 30, 30, 
        31, 32, 32, 33, 33, 33, 33, 32, 33, 33, 33, 33, 33, 33, 33, 33, 36, 34, 35, 
        34, 35, 38, 37, 38, 38, 39, 38, 41, 39, 43, 38, 41, 42, 43, 43, 40, 41, 42, 
        43, 44, 45, 46, 47, 48, 49, 50, 48, 49, 53, 51, 51, 51, 53, 55, 56, 55, 56, 
        58, 58, 58, 59, 61, 61, 63, 63, 64, 64, 64, 65, 65, 65, 64, 64, 65, 65, 65, 
        66, 67, 65, 66, 65, 68, 66, 65, 66, 65, 66, 67, 65, 66, 67, 68, 69, 70, 71, 
        72, 71, 72, 71, 76, 71, 76, 72, 71, 76, 71, 73, 72, 76, 76, 73, 73, 72, 76, 
        76, 75, 76, 76, 75, 81, 81, 83, 81 };

    splitsB = new int[] { 
        0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 
        12, 12, 12, 12, 12, 12, 12, 13, 14, 14, 15, 16, 16, 16, 16, 17, 16, 19, 17, 
        19, 19, 19, 20, 21, 22, 22, 22, 22, 22, 22, 22, 24, 22, 22, 24, 24, 26, 27, 
        27, 28, 28, 28, 30, 28, 28, 28, 29, 28, 28, 28, 29, 29, 30, 33, 33, 33, 33, 
        34, 35, 37, 37, 37, 38, 38, 38, 37, 38, 38, 38, 38, 38, 39, 38, 44, 43, 44, 
        39, 40, 41, 43, 43, 43, 45, 46, 46, 46, 47, 48, 48, 49, 49, 50, 51, 51, 49, 
        49, 50, 51, 50, 51, 50, 50, 51, 50, 51, 51, 51, 53, 55, 55, 55, 56, 56, 56, 
        56, 56, 57, 58, 61, 61, 61, 63, 63, 63, 64, 65, 66, 65, 66, 66, 66, 65, 66, 
        66, 66, 66, 66, 68, 71, 66, 66, 68, 67, 71, 66, 66, 68, 67, 71, 66, 66, 68, 
        68, 71, 70, 70, 72, 72, 76, 75, 75, 78, 78, 78, 80, 80, 80, 80, 81, 81, 81, 
        82, 83, 84, 85, 86, 86, 86, 86, 86, 86, 88, 86, 90, 86, 92, 87, 86, 89, 86, 
        92, 87, 86, 87, 86, 91, 89, 89, 90, 90, 92, 92, 92, 93, 93, 93, 95, 95, 95, 
        95, 95, 95, 95, 95 };

    baseCurve = new MontCurve();
  }


  public int geteB() {
    return eB;
  }


  public int geteA() {
    return eA;
  }


  public BigInteger getOrderA() {
    return orderA;
  }


  public BigInteger getOrderB() {
    return orderB;
  }


  public MontCurve getCurve() {
    return new MontCurve (baseCurve);
  }


  public FpPointAffine getGeneratorA() {
    return pA;
  }


  public FpPointAffine getGeneratorB() {
    return pB;
  }


  public int getMIPA() {
    return maxIntPointsA;
  }


  public int getMIPB() {
    return maxIntPointsB;
  }


  public int getMaxA() {
    return maxA;
  }


  public int getMaxB() {
    return maxB;
  }


  public int[] getSplitsA() {
    int[] retval = new int[maxA];
    System.arraycopy (splitsA, 0, retval, 0, maxA);
    return retval;
  }


  public int[] getSplitsB() {
    int[] retval = new int[maxB];
    System.arraycopy (splitsB, 0, retval, 0, maxB);
    return retval;
  }


  public int getObits() {
    return obits;
  }


  public SidhKeyPair generateKeyPair (int aOrB) {
    // Generate a key pair with a randomly generated private key
    SidhKeyPair keys = new SidhKeyPair (aOrB, this);
    return keys;
  }


  public SidhKeyPair generateKeyPair (int aOrB, BigInteger prKey) {
    // Generate a key pair with the specified private key
    SidhKeyPair keys = new SidhKeyPair (aOrB, prKey, this);
    return keys;
  }


  public byte[] calculateAgreementA (SidhPrivateKey privKeyA, SidhPublicKey pubKeyB) {
    // A's shared secret generation
    // Inputs: A's private key, B's public key
    // Outputs: a shared secret which is an element of the quadratic extension field Fp^2

    int i, row, m, index = 0, ptsIdx[], npts = 0;
    F2Point r, pts[];

    // Extract B's public key
    FourIsogeny fourIsog = new FourIsogeny (pubKeyB.getA());
    F2elm pkB2 = pubKeyB.getP();
    F2elm pkB3 = pubKeyB.getQ();
    F2elm pkB4 = pubKeyB.getD();

    pts = new F2Point[maxIntPointsA];
    ptsIdx = new int[maxIntPointsA];

    r = fourIsog.ladder3pt (pkB2, pkB3, pkB4, privKeyA.getKey(), obits);
    r = fourIsog.first4Isog (r);

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

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    fourIsog.get4Isog (r);
    return fourIsog.jInv().toByteArray();
  }


  public byte[] calculateAgreementB (SidhPrivateKey privKeyB, SidhPublicKey pubKeyA) {
    // B's shared secret generation
    // Inputs: B's private key, A's public key
    // Outputs: a shared secret which is an element of the quadratic extension field Fp^2

    int i, row, m, index = 0, ptsIdx[], npts = 0;
    F2Point r, pts[];
    F2elm pkA2, pkA3, pkA4;
    ThreeIsogeny threeIsog;

    pts = new F2Point[maxIntPointsB];
    ptsIdx = new int[maxIntPointsB];

    // Extract A's public key
    threeIsog = new ThreeIsogeny (pubKeyA.getA());
    pkA2 = pubKeyA.getP();
    pkA3 = pubKeyA.getQ();
    pkA4 = pubKeyA.getD();

    r = threeIsog.ladder3pt(pkA2, pkA3, pkA4, privKeyB.getKey(), obits);

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

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    threeIsog.get3Isog (r);
    return threeIsog.jInv().toByteArray();
  }
}
