
package sidh;

/**************************************************************************************************
 *
 * Implements private keys for the Supersingular Isogeny Diffie-Hellman key exchange algorithm. 
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

class SidhPrivateKey {
  private final BigInteger m;           


  public SidhPrivateKey (BigInteger mIn) {
    m = mIn;
  }


  public SidhPrivateKey (SidhPrivateKey kIn) {
    m = kIn.getKey ();
  }


  public SidhPrivateKey (byte[] bytesIn) {
    m = new BigInteger (bytesIn);
  }


  public SidhPrivateKey (int aOrB, BigInteger order) {
    // Generate a random private key
    BigInteger temp, randmod, three = BigInteger.valueOf(3);
    boolean condition;

    temp = Felm.genRandom (order);
    if (aOrB == SidhKeyExchange.PARTYA) 
      condition = temp.testBit(0);
    else { 
      randmod = temp.mod(three);
      condition = randmod.intValue() != 0;  
    }

    while (temp.equals(BigInteger.ZERO) || condition) {
        temp = Felm.genRandom (order);
        if (aOrB == SidhKeyExchange.PARTYA) 
          condition = temp.testBit(0);
        else { 
          randmod = temp.mod(three);
          condition = randmod.intValue() != 0;  
        }
    }

    m = temp;
  }


  public BigInteger getKey() {
    return m;
  }    


  public boolean privateKeyEquals (SidhPrivateKey otherKey) {
    return m.compareTo (otherKey.getKey ()) == 0;
  }
    

  public byte[] serialize() {
    return m.toByteArray();
  }
}
