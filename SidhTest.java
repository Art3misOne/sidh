
package sidh;

/**************************************************************************************************
 *
 * Basic tests for supersingular isogeny Diffie-Hellman key exchange and underlying field and curve
 * arithmetic.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;

class SidhTest {
  public static boolean testrandom = true;

  public static void main (String[] args) {
    // Using default parameters
    SidhKeyExchange kex = new SidhKeyExchange();

    SidhKeyPair keysA, keysB;
    SidhPublicKey tempkey;
    byte[] sharedA, sharedB;

    if (testrandom) {  
      System.out.println ("\nTesting key exchange/validation with randomly generated private keys\n");

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
