package rlwe;

/**************************************************************************************************
 *
 * Basic tests for the Ring Learning With Errors key exchange.
 *
 **************************************************************************************************/

import java.math.BigInteger;


class RLWETest {

  public static void main (String[] args) {
    /**** Tests for Fq arithmetic ****/
    Felm a, b;

    Felm.setModulus (2147483631);
    System.out.println ("Modulus set to " + Felm.getModulus()); 

    a = new Felm ();
    b = new Felm (-7);

    a.fqSetValue (253);
    System.out.println ("a set to " + a.fqGetValue());
    System.out.println ("b set to " + b.fqGetValue() + "\n");

    System.out.println ("a + b = " + a.fqAdd (b));
    System.out.println ("b - a = " + b.fqSub (a));
    System.out.println ("a * b = " + a.fqMult (b));
    System.out.println ("a**2 = " + a.fqSqr() + "\n");

    System.out.println ("a is zero? " + a.fqIsZero());
    System.out.println ("ZERO is zero? " + Felm.ZERO.fqIsZero());
    System.out.println ("a is even? " + a.fqIsEven());
    System.out.println ("b is even? " + b.fqIsEven());
    System.out.println ("a == b? " + a.fqEqual(b)); 
    System.out.println ("a < b? " + a.fqLessThan(b));
    System.out.println ("a > b? " + a.fqGreaterThan(b) + "\n");

    System.out.println ("-a = " + a.fqNegate());
    System.out.println ("a**-1 = " + a.fqInverse());
    System.out.println ("a / 2 = " + a.fqDiv2());
    System.out.println ("a / b = " + a.fqDiv(b) + "\n");

    System.out.println ("a as a string = " + a.toString());
    System.out.println ("a as a byte array = " + a.toByteArray() + "\n");

    /**** Tests for Polynomial ****/
    long[] fvalues = new long [] {3, 16, 1, 8, 17};
    long[] gvalues = new long [] {6, 8, 13, 2, 18};
    long[] hvalues = new long [] {1, 1, 1, 1}; 
    Polynomial f = new Polynomial (fvalues);
    Polynomial g = new Polynomial (gvalues);
    Polynomial h = new Polynomial (hvalues);

    System.out.println ("f = " + f.toString());
    System.out.println ("g = " + g);
    System.out.println ("h = " + h + "\n");
    
    System.out.println ("Degree of f = " + f.getDegree());
    System.out.println ("f[3] = " + f.getCoeff(3));
    f.setCoeff (3, 20);
    System.out.println ("Changed f[3] to " + f.getCoeff(3));
    System.out.println ("f = " + f + "\n");

    System.out.println ("f == f? " + f.equal(f));
    System.out.println ("f == g? " + f.equal(g) + "\n");

    System.out.println ("f + h = " + f.polyAdd(h));
    System.out.println ("f - g = " + f.polySub(g));
    System.out.println ("f * g = " + f.polyMult(g));
    System.out.println ("f*g + h = " + f.polyMultAdd (g, h));
    System.out.println ("f % h = " + f.polyMod (h));

    System.out.println ("f as a byte array: " + f.toByteArray() + "\n"); 
    

    /**** Tests for elements of a polynomial quotient ring ****/
    RingElt.setModulus (h);
    RingElt rf = new RingElt (f);
    RingElt rg = new RingElt (g);
    
    System.out.println (rf);
    System.out.println (rg);

    System.out.println ("rf + rg = " + rf.ringAdd(rg));
    System.out.println ("rf - rg = " + rf.ringSub(rg));
    System.out.println ("rf * rg = " + rf.ringMult(rg) + "\n");


    /**** Sampling timing tests ****/

    int iterations = 1000, i;
    long startTime, endTime, gavg = 0, gctavg = 0, bavg = 0;
    long[] gtimes, gcttimes, btimes;

    gtimes = new long[iterations];
    gcttimes = new long[iterations];
    btimes = new long[iterations];

    for (i = 0; i < iterations; i++) {
      startTime = System.nanoTime();
      rf = RLWEPrivateKey.sampleGaussian(1024);
      endTime = System.nanoTime();

      gtimes[i] = (endTime - startTime) / 1000;        // Convert nanoseconds to microseconds
      gavg += gtimes[i];

      startTime = System.nanoTime();
      rg = RLWEPrivateKey.sampleGaussianCT(1024);    
      endTime = System.nanoTime();

      gcttimes[i] = (endTime - startTime) / 1000;      // Convert nanoseconds to microseconds
      gctavg += gcttimes[i];

      startTime = System.nanoTime();
      rg = RLWEPrivateKey.sampleBinomial(1024);    
      endTime = System.nanoTime();

      btimes[i] = (endTime - startTime) / 1000;        // Convert nanoseconds to microseconds
      bavg += btimes[i];
    }

    gavg = gavg / iterations;
    gctavg = gctavg / iterations;
    bavg = bavg / iterations;

    long gvariance = 0, gctvariance = 0, bvariance = 0;

    for (i = 0; i < iterations; i++) {
	gvariance += Math.pow (gtimes[i] - gavg, 2); 
        gctvariance += Math.pow (gcttimes[i] - gctavg, 2);
	bvariance += Math.pow (btimes[i] - bavg, 2); 
    }

    gvariance = gvariance / iterations;
    gctvariance = gctvariance / iterations;
    bvariance = bvariance / iterations;

    System.out.println ("With times given in microseconds: ");
    System.out.print ("Average time for one Gaussian sample: " + gavg + " with variance = ");
    System.out.println (gvariance + " (std dev = " + Math.sqrt(gvariance) + ")"); 
    System.out.print ("Average constant time Gaussian sample: " + gctavg + " with variance = ");
    System.out.println (gctvariance + " (std dev = " + Math.sqrt(gctvariance) + ")"); 
    System.out.print ("Average time for one Binomial sample: " + bavg + " with variance = ");
    System.out.println (bvariance + " (std dev = " + Math.sqrt(bvariance) + ")\n"); 


    /**** Tests for keys ****/

    RLWEKeyPair keysI, keysR;

    keysI = new RLWEKeyPair (rf);
    keysR = new RLWEKeyPair (rf);

    System.out.println ("Testing key generation with toy sized keys\n");

    System.out.println ("Initiator private key: ");
    System.out.println ("\t s = " + keysI.getPrivateKey().getS());
    System.out.println ("\t e = " + keysI.getPrivateKey().getE() + "\n");

    System.out.println ("Initiator public key: ");
    System.out.println ("\t k = " + keysI.getPublicKey().getKey() + "\n");

    System.out.println ("Responder private key: ");
    System.out.println ("\t s = " + keysR.getPrivateKey().getS());
    System.out.println ("\t e = " + keysR.getPrivateKey().getE() + "\n");

    System.out.println ("Responder public key: ");
    System.out.println ("\t k = " + keysR.getPublicKey().getKey() + "\n");
    

    /**** Tests for key exchange ****/

    RLWEKeyExchange kex = new RLWEKeyExchange ();
    BigInteger recData, secretI, secretR;
    BigInteger[] response = new BigInteger[2];

    
    keysI = kex.generateKeyPair ();
    keysR = kex.generateKeyPair ();

    System.out.println ("Testing key exchange with real sized keys\n");

    /*
    System.out.println ("Initiator private key: ");
    System.out.println ("\t s = " + keysI.getPrivateKey().getS());
    System.out.println ("\t e = " + keysI.getPrivateKey().getE() + "\n");

    System.out.println ("Initiator public key: ");
    System.out.println ("\t k = " + keysI.getPublicKey().getKey() + "\n");

    System.out.println ("Responder private key: ");
    System.out.println ("\t s = " + keysR.getPrivateKey().getS());
    System.out.println ("\t e = " + keysR.getPrivateKey().getE() + "\n");

    System.out.println ("Responder public key: ");
    System.out.println ("\t k = " + keysR.getPublicKey().getKey() + "\n");
    */

    response = kex.respAgreement (keysI.getPrivateKey(), keysR.getPublicKey());
    recData = response[0];
    secretR = response[1];

    System.out.println ("Reconciliation data = " + recData + "\n");

    secretI = kex.initAgreement (keysR.getPrivateKey(), keysI.getPublicKey(), recData);

    if (secretI.equals(secretR)) {
      System.out.println ("Shared secrets match :)\n");
      String debug = "Shared secret = " + secretI + "\n";
      System.out.println (debug);
    }

    else {
      System.out.println ("Shared secrets do not match :(\n");
      String debug = "secretI = " + secretI + "\nsecretR = " + secretR + "\n";
      System.out.println (debug);
    }
  }
}
