
package rlwe;

/**************************************************************************************************
 *
 * Implements objects for public and private keys including key generation functions.
 *  
 * TODO: 
 *    - write sampling functions from various error distributions
 *
 **************************************************************************************************/

import java.util.Random;
import java.math.BigInteger;
import java.util.Arrays;

class RLWEPublicKey {
  RingElt key;


  public RLWEPublicKey (RLWEPrivateKey k, RingElt a) {
    key = a.ringMultAdd(k.getS (), k.getE ());
  }


  public RLWEPublicKey (RingElt b) {
    key = new RingElt (b);
  }

  
  public RLWEPublicKey (RLWEPublicKey k) {
    key = new RingElt (k.key);
  }


  public RLWEPublicKey (byte[] inBytes) {
    key = new RingElt (inBytes);
  }


  public RingElt getKey () {
    return key;
  }


  public byte[] serialize () {
    return key.toByteArray();
  }


  public boolean equal (RLWEPublicKey pk) {
    return key.equal (pk.key);
  }


  public int hashcode () {
    return Arrays.hashCode (serialize());
  }
}


class RLWEPrivateKey {
  private final RingElt s;
  private final RingElt e;
  private final RingElt eprime;


  public RLWEPrivateKey (RingElt sIn, RingElt eIn) {
    s = new RingElt (sIn);
    e = new RingElt (eIn);
    eprime = sampleBinomial (s.getDegree());
  }


  public RLWEPrivateKey (int degree) {
    s = sampleBinomial (degree);
    e = sampleBinomial (degree);
    eprime = sampleBinomial (degree);
  }


  public RLWEPrivateKey (byte[] inBytes) {
    // Reconstruct a private key from a byte array assuming s and e are the same size.
    int len = inBytes.length;
    s = new RingElt (Arrays.copyOfRange (inBytes, 0, len / 2)); 
    e = new RingElt (Arrays.copyOfRange (inBytes, len / 2, len));
    eprime = null;
  }


  public RingElt getS () {
    return s;
  }
  

  public RingElt getE () {
    return e;
  }


  public RingElt getEprime () {
    return eprime;
  }


  public byte[] serialize () {
    byte[] sba = s.toByteArray ();
    byte[] eba = e.toByteArray ();
    byte[] r = new byte[sba.length + eba.length];
    System.arraycopy (sba, 0, r, 0, sba.length);
    System.arraycopy (eba, 0, r, sba.length, eba.length);
    return r;
  }


  public static long singleSample (BigInteger in) {
    int i = 0;

    while (in.compareTo (RLWEUtilities.GAUSSIAN_TABLE[i]) >= 0)
      i++;

    return (long) i;
  }


  public static long selectCT (long x, long y, boolean selector) {
    long mask = 0 - Boolean.compare (selector, false);   // use compare to convert boolean to long
    return (x & mask) | (y & ~mask);
  }


  public static boolean lessThanCT (BigInteger x, BigInteger y) {
    return (x.subtract(y)).signum() == -1; 
  }


  public static long singleSampleCT (BigInteger in) {
    long index = 0;
    int tablelen = RLWEUtilities.GAUSSIAN_TABLE.length;

    for (int i = 0; i < tablelen; i++)
      index = selectCT (index, (long) i+1, lessThanCT (RLWEUtilities.GAUSSIAN_TABLE[i], in));

    return index;
  }


  public static RingElt sampleGaussian (int degree) {
    int i, m;
    long[] s = new long[degree];
    Random rand = new Random ();
    BigInteger rnd, r = new BigInteger (degree, rand);

    for (i = 0; i < degree; i++) {
      rnd = new BigInteger (192, rand);
      m = boolToInt (r.testBit (i));
      s[i] = singleSample (rnd); 
      if (m != 0) 
        s[i] = -s[i];
    }

    return new RingElt (s);
  }


  public static boolean isNonzeroCT (long x) {
    return ((x | -x) >> (Long.SIZE - 1)) == 1;
  }


  public static RingElt sampleGaussianCT (int degree) {
    int i, m;
    long[] s = new long[degree];
    Random rand = new Random ();
    BigInteger rnd, r = new BigInteger (degree, rand);

    for (i = 0; i < degree; i++) {
      rnd = new BigInteger (192, rand);
      m = boolToInt (r.testBit (i));
      s[i] = singleSampleCT (rnd); 
      s[i] = selectCT (s[i], -s[i], isNonzeroCT (m));
    }

    return new RingElt (s);
  }  


  public static RingElt sampleBinomial (int degree) {
    int i, j, idx, b0, b1, iters = RLWEUtilities.BINOMITERATIONS;
    long[] s = new long[degree];
    Random rand = new Random ();
    BigInteger randbits0 = new BigInteger (degree*iters, rand);
    BigInteger randbits1 = new BigInteger (degree*iters, rand);

    for (i = 0; i < degree; i++) {
      for (j = 0; j < iters; j++) {
        idx = i * iters + j;
        b0 = boolToInt (randbits0.testBit (idx));
        b1 = boolToInt (randbits1.testBit (idx));
	s[i] += b1 - b0;
      }
    }

    return new RingElt (s);
  }  


  public static int boolToInt (boolean b) {
    if (b) return 1;
    return 0;
  }
}


class RLWEKeyPair {
  private final RLWEPublicKey pubKey; 
  private final RLWEPrivateKey privKey;


  public RLWEKeyPair (RLWEPrivateKey prKey, RingElt a) {
    privKey = prKey;
    pubKey = new RLWEPublicKey (prKey, a);
  }


  public RLWEKeyPair (RingElt a) {
    privKey = new RLWEPrivateKey (a.degree);
    pubKey = new RLWEPublicKey (privKey, a);
  }


  public RLWEPrivateKey getPrivateKey () {
    return privKey;
  }


  public RLWEPublicKey getPublicKey () {
    return pubKey;
  }
}
