
package sidh;

/**************************************************************************************************
 *
 * Implements public/private keypair for the Supersingular Isogeny Diffie-Hellman key exchange 
 * algorithm. 
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;

class SidhKeyPair {
  private final SidhPublicKey pubKey;
  private final SidhPrivateKey privKey;


  public SidhKeyPair (int aOrB, BigInteger prKey, SidhKeyExchange kex) {
    // Create a key-pair where the private key = prKey
    privKey = new SidhPrivateKey (prKey);
    pubKey = new SidhPublicKey (aOrB, privKey, kex);
  }


  public SidhKeyPair (int aOrB, SidhKeyExchange kex) {
    // Create a key-pair where the private key is generated randomly
    BigInteger order;

    if (aOrB == SidhKeyExchange.PARTYA)
      order = kex.getOrderA();
    else
      order = kex.getOrderB();

    privKey = new SidhPrivateKey (aOrB, order);
    pubKey = new SidhPublicKey (aOrB, privKey, kex);
  }


  public SidhKeyPair (SidhPublicKey publicK, SidhPrivateKey privateK) {
    pubKey = publicK;
    privKey = privateK;
  }


  public SidhPublicKey getPublicKey() {
    return pubKey;
  }

    
  public SidhPrivateKey getPrivateKey() {
    return privKey;
  } 
}
