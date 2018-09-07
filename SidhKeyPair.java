
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


  public SidhKeyPair (int aOrB, SidhPrivateKey prKey, SidhKeyExchange kex) {
    privKey = new SidhPrivateKey (prKey);
    pubKey = new SidhPublicKey (aOrB, privKey, kex);
  }


  public SidhKeyPair (int aOrB, SidhKeyExchange kex) {
    BigInteger order;

    if (aOrB == SidhKeyExchange.PARTYA)
      order = kex.getOrderA ();
    else
      order = kex.getOrderB ();
    
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
