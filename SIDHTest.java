
package sidh;

/**************************************************************************************************
 *
 * Basic tests for supersingular isogeny Diffie-Hellman key exchange and underlying field and curve
 * arithmetic.
 *  
 **************************************************************************************************/

import java.math.BigInteger;

class SIDHTest {
  public static boolean testrandom = true;

  public static void main (String[] args) {
    // Using default parameters
    SIDHKeyExchange kex = new SIDHKeyExchange();

    SIDHKeyPair keysA, keysB;
    SIDHPublicKey tempkey;
    F2elm sharedA, sharedB;

    if (testrandom) {  
      System.out.println ("\nTesting key exchange/validation with randomly generated private keys\n");

      keysA = kex.generateKeyPair (SIDHKeyExchange.PARTYA);

      System.out.println ("A's Keys:");
      System.out.println ("\tPrivate: " + keysA.getPrivateKey().getKey());
      tempkey = keysA.getPublicKey();
      System.out.println ("\tPublic: " + tempkey.getA());
      System.out.println ("\t        " + tempkey.getP());
      System.out.println ("\t        " + tempkey.getQ());
      System.out.println ("\t        " + tempkey.getD() + "\n");

      int eA = kex.geteA(), eB = kex.geteB();

      SIDHKeyValidate validator = new SIDHKeyValidate (keysA.getPublicKey());
      if (validator.validateAsKey (eA, eB))
        System.out.println ("A's public key is valid\n");
      else 
        System.out.println ("A's public key is not valid :(\n");

      keysB = kex.generateKeyPair (SIDHKeyExchange.PARTYB);

      System.out.println ("B's Keys (original, then serialized and reconstructed):");
      System.out.println ("\tPrivate: " + keysB.getPrivateKey().getKey());
      SIDHPrivateKey reconstruct = new SIDHPrivateKey (keysB.getPrivateKey().serialize());
      System.out.println ("\tReconst: " + reconstruct.getKey() + "\n");
      tempkey = keysB.getPublicKey();
      System.out.println ("\tPublic: " + tempkey.getA());
      System.out.println ("\t        " + tempkey.getP());
      System.out.println ("\t        " + tempkey.getQ());
      System.out.println ("\t        " + tempkey.getD() + "\n");
      tempkey = new SIDHPublicKey (tempkey.serialize());
      System.out.println ("\tRecon:  " + tempkey.getA());
      System.out.println ("\t        " + tempkey.getP());
      System.out.println ("\t        " + tempkey.getQ());
      System.out.println ("\t        " + tempkey.getD());

      validator = new SIDHKeyValidate (keysB.getPublicKey());
      if (validator.validateBsKey (eA, eB))
        System.out.println ("B's public key is valid\n");
      else
	System.out.println ("B's public key is not valid\n");

      sharedA = kex.calculateAgreementA (keysA.getPrivateKey(), keysB.getPublicKey());
      sharedB = kex.calculateAgreementB (keysB.getPrivateKey(), keysA.getPublicKey());

      if (sharedA.f2Equals (sharedB)) {
        System.out.println ("Shared secrets match :)\n");
        String debug = "Shared secret = " + sharedA + "\n";
        System.out.println (debug);
      }

      else {
        System.out.println ("Shared secrets do not match :(\n");
        String debug = "SharedA = " + sharedA + "\nSharedB = " + sharedB;
        System.out.println (debug + "\n");
      }
    }
  }
}
