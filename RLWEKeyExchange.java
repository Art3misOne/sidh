package rlwe;

/**************************************************************************************************
 *
 * Implements RLWE key exchange algorithm.
 *
 * TODO:
 *    - Create multiple constructors based on parameter sets
 *    - Implement reconciliation function and function to compute reconciliation data
 * 
 **************************************************************************************************/

import java.util.Arrays;
import java.math.BigInteger;
import java.util.Random;

class RLWEKeyExchange {
  public int m;                       // degree of the polynomial
  public long q;                      // modulus
  public RingElt a;                   // parameter for a*s + e computations
  public byte recMethod;              // reconciliation method


    public RLWEKeyExchange (int mIn, long qIn, byte rm) {
    m = mIn;
    q = qIn;
    Felm.setModulus (q);
    setQuotCyclotomic ();
    a = RLWEUtilities.setA (m, q);
    recMethod = rm;
  }


  public RLWEKeyExchange () {
    // Use a default set of parameters if none is provided
    m = 1024;
    q = 40961;
    Felm.setModulus (q);
    setQuotCyclotomic ();
    a = RLWEUtilities.setA (m, q);
    recMethod = RLWEUtilities.PEIKERT;
  }


  public RingElt getA () {
    return a;
  }


  public void setQuotCyclotomic () {
    // Set the modulus of the quotient ring to be the cyclotomic polynomial
    // Assumes m is either an odd prime or a power of 2
    long[] qucoef = new long[m];

    if (m%1 == 1) 
      Arrays.fill (qucoef, 1); 
    else
      qucoef[0] = qucoef[m-1] = 1;

    RingElt.setModulus (new Polynomial (qucoef));
  }


  public void setQuotient (Polynomial qq) {
    RingElt.setModulus (qq);
  }


  public RLWEKeyPair generateKeyPair () {
    return new RLWEKeyPair (a);
  }


  public BigInteger[] roundCrossRound (RingElt k) {
    // This function does simultaneous round and cross round. Before doing so, it performs a 
    // random nudge for the edge cases to remove bias.

    int i, d = k.getDegree();
    Random rnd = new Random ();
    BigInteger r = new BigInteger (d, rnd);      // Get one bit of random per coefficient
    BigInteger round = BigInteger.ZERO, crossround = BigInteger.ZERO;
    Felm edge, coeff;
    long edgeLong, cLong;

    // Set edge case
    edgeLong = q/4;
    if (q % 4 != 1)
      edgeLong *= 3;
    edge = new Felm (edgeLong - 1);

    for (i = 0; i < d; i++) {
      coeff = k.getCoeff (i);

      if (r.testBit (i))
        coeff = randomNudge (coeff, edge);
      r = r.shiftRight (1);

      cLong = coeff.fqGetValue ();

      if ((q/4) < cLong && cLong < (q*3/4))
        round = round.setBit (i);
      if ( ((q/4) < cLong && cLong <= (q/2)) || (q*3/4) <= cLong )
        crossround = crossround.setBit(i);
    }

    return new BigInteger[] { crossround, round };
  }


  private Felm randomNudge (Felm coeff, Felm edge) {
    if (coeff.fqEqual (Felm.ZERO))
      return coeff.fqSub (Felm.ONE);
    if (coeff.fqEqual (edge))
      return edge.fqAdd (Felm.ONE);
    return coeff;
  }


  public BigInteger recPeikert (RingElt k, BigInteger recData) {
    double lower0 = Math.round ((3.0/8.0) * q), lower1 = Math.round ((1.0/8.0) * q);
    double upper0 = Math.round ((7.0/8.0) * q), upper1 = Math.round ((5.0/8.0) * q);
    double cdbl;
    int i, d = k.getDegree ();
    BigInteger key = BigInteger.ZERO;

    for (i = 0; i < d; i++) {
      cdbl = (double) k.getCoeff (i).fqGetValue ();

      if (recData.testBit (i)) {
        if (lower1 < cdbl && cdbl < upper1)
          key = key.setBit (i);
      }

      else {
        if (lower0 < cdbl && cdbl < upper0)
          key = key.setBit (i);   
      }
    }

    return key;
  }


  public BigInteger initAgreement (RLWEPrivateKey kI, RLWEPublicKey kR, BigInteger w) {
    // Input: Initiator's priv key, Responder's pub key, reconciliation data, reconciliation method
    // Output: shared secret

    RingElt k = kR.getKey().ringMultAdd (kI.getS(), kI.getEprime());
    if (recMethod == RLWEUtilities.PEIKERT) 
      return recPeikert (k, w);
    else 
      return BigInteger.ZERO;
  }


  public BigInteger[] respAgreement (RLWEPrivateKey kR, RLWEPublicKey kI) {
    // Input: Responder's private key, Initiator's public key, reconciliation method
    // Output: { shared secret, reconciliation data }

    RingElt k = kI.getKey().ringMultAdd (kR.getS(), kR.getEprime());
    if (recMethod == RLWEUtilities.PEIKERT)
      return roundCrossRound (k);
    else
      return new BigInteger[] { BigInteger.ZERO, BigInteger.ZERO };
  }
}
