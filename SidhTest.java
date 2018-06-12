
package sidh;

/**************************************************************************************************
 *
 * Basic tests for supersingular isogeny Diffie-Hellman key exchange and validation.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;
import java.security.SecureRandom;


class SidhTest {
  public static boolean testrandom = true;

  public static void main (String[] args) {
    // Arithmetic tests
    F2elm a, asq, asqsqrt, asq2;
    Felm a0;
    BigInteger p, dandr[];
    SecureRandom rnd = new SecureRandom();

    for (int i = 0; i < 100; i++) {
      p = new BigInteger (100, 10, rnd);
	
      try {
	Felm.setPrime(p);
      } catch (InvalidFieldException e) {
	System.exit(-1);
      }

      a = new F2elm (rnd);

      asq = a.f2Sqr ();                // asq constructed to be a quadratic residue
      asqsqrt = asq.f2Sqrt ();
      asq2 = asqsqrt.f2Sqr ();         // sqrt (a^2) might yield a diff residue than a

      if (asq2.f2Equals (asq) == false) {
        System.out.println ("Square root test iteration " + i + ": ");
	System.out.println ("a = " + a + "\t\t asq = " + asq);
	System.out.println ("asqsqrt = " + asqsqrt + "\t asq2 = " + asq2 + "\n");
	break;
      }
    }
    
    // Using default parameters
    SidhKeyExchange kex = new SidhKeyExchange();          
    SidhKeyPair keysA, keysB;
    SidhPublicKey tempkey;
    byte[] sharedA, sharedB;

    Felm fortyseven = new Felm (47);
    Felm fiftytwo = new Felm (52);
    Felm four = new Felm (4);


    if (testrandom) {  
      System.out.println ("\nTesting exchange/validation with randomly generated private keys\n");

      keysA = kex.generateKeyPair (SidhKeyExchange.PARTYA);

      //System.out.println ("A's Keys:");
      //System.out.println ("\tPrivate: " + keysA.getPrivateKey().getKey());
      //tempkey = keysA.getPublicKey();
      //System.out.println ("\tPublic: " + tempkey.getA());
      //System.out.println ("\t        " + tempkey.getP());
      //System.out.println ("\t        " + tempkey.getQ());
      //System.out.println ("\t        " + tempkey.getD() + "\n");

      int eA = kex.geteA(), eB = kex.geteB();

      SidhKeyValidate validator = new SidhKeyValidate (keysA.getPublicKey());
      if (validator.validateAsKey (eA, eB))
        System.out.println ("\nA's public key is valid");
      else 
        System.out.println ("\nA's public key is not valid :(");

      keysB = kex.generateKeyPair (SidhKeyExchange.PARTYB);

      //System.out.println ("B's Keys (original, then serialized and reconstructed):");
      //System.out.println ("\tPrivate: " + keysB.getPrivateKey().getKey());
      //SidhPrivateKey reconstruct = new SidhPrivateKey (keysB.getPrivateKey().serialize());
      //System.out.println ("\tReconst: " + reconstruct.getKey() + "\n");
      //tempkey = keysB.getPublicKey();
      //System.out.println ("\tPublic: " + tempkey.getA());
      //System.out.println ("\t        " + tempkey.getP());
      //System.out.println ("\t        " + tempkey.getQ());
      //System.out.println ("\t        " + tempkey.getD() + "\n");
      //tempkey = new SidhPublicKey (tempkey.serialize());
      //System.out.println ("\tRecon:  " + tempkey.getA());
      //System.out.println ("\t        " + tempkey.getP());
      //System.out.println ("\t        " + tempkey.getQ());
      //System.out.println ("\t        " + tempkey.getD());

      validator = new SidhKeyValidate (keysB.getPublicKey());
      if (validator.validateBsKey (eA, eB))
        System.out.println ("B's public key is valid\n");
      else
	System.out.println ("B's public key is not valid\n");

      sharedA = kex.calculateAgreementA (keysA.getPrivateKey(), keysB.getPublicKey());
      sharedB = kex.calculateAgreementB (keysB.getPrivateKey(), keysA.getPublicKey());

      if (Arrays.equals (sharedA, sharedB)) {
        System.out.print ("Shared secrets match :)\n\n Shared secret = ");
        printByteArray (sharedA);
        System.out.println ("\n");
      }

      else {
        System.out.print ("Shared secrets do not match :(\n\nSharedA = ");
        printByteArray (sharedA);
        System.out.print ("\nSharedB = ");
        printByteArray (sharedB);
        System.out.println ("\n");
      }
    }
  }

  public static void printByteArray (byte[] in) {
    System.out.print ("0x");
    for (int i = 0; i < in.length; i++)
      System.out.printf ("%02x ", in[i]);
  }
}
