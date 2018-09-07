
package sidh;

/**************************************************************************************************
 *
 * Basic tests for supersingular isogeny Diffie-Hellman key exchange and validation.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;
import java.security.SecureRandom;
import java.lang.management.*;

class SidhTest {
  public static boolean testfixed = true;

  public static void main (String[] args) {
    // Using default parameters
    SidhKeyExchange kex = new SidhKeyExchange("sidhP503");
    SidhKeyPair keysA, keysB;
    SidhPublicKey reconstructedApub, publicB;
    SidhPrivateKey reconstructedApriv, privateB;
    byte[] sharedA, sharedB, pubKeyBytes, privKeyBytes;
    
    long startTime, startTime2, endTime, totalTime = 0;
    int i, iterations = 10;
    BigInteger akey, bkey;
    
    // Testing key exchange
    
    if (testfixed) {
      akey = new BigInteger ("2b701ec1698bf9a513875fb7188c1d63fbd59ac8a378c3fbb1c98496173f6e", 16);
      bkey = new BigInteger ("9cfe2a283dfb23c330fb2202dd2c34f8a0c45f2ab761ec7ca4bc11a3324d5c7", 16);
      
      keysA = kex.generateKeyPair (SidhKeyExchange.PARTYA, new SidhPrivateKey (akey));
      keysB = kex.generateKeyPair (SidhKeyExchange.PARTYB, new SidhPrivateKey (bkey));

      pubKeyBytes = keysA.getPublicKey().serialize();
      privKeyBytes = keysA.getPrivateKey().serialize();

      reconstructedApub = new SidhPublicKey (pubKeyBytes);
      reconstructedApriv = new SidhPrivateKey (privKeyBytes);

      if (reconstructedApub.publicKeyEquals (keysA.getPublicKey()))
        System.out.println ("Public key reconstruction successful");
      else
        System.out.println ("Public key reconstruction unsuccessful");

      if (reconstructedApriv.privateKeyEquals (keysA.getPrivateKey()))
        System.out.println ("Private key reconstruction successful\n");
      else
        System.out.println ("Private key reconstruction unsuccessful\n");      
      
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

    if (iterations > 0) {
      long mintime = 10000000;
      long maxtime = 0;
      long thistime;
      
      System.out.println ("\nRunning timing tests\n");

      keysB = kex.generateKeyPair (SidhKeyExchange.PARTYB);
      publicB = keysB.getPublicKey ();
      privateB = keysB.getPrivateKey ();
      
      for (i = 0; i < iterations; i++) {
	startTime = getCpuTime ();
	  
	keysA = kex.generateKeyPair (SidhKeyExchange.PARTYA);
	sharedA = kex.calculateAgreementA (keysA.getPrivateKey(), publicB);
	sharedB = kex.calculateAgreementB (privateB, keysA.getPublicKey());

	if (Arrays.equals (sharedA, sharedB) == false) {
	    System.out.println ("Failed exchange at iteration " + i + "\n");
	    break;
	}
	
	endTime = getCpuTime();
	thistime = (endTime - startTime) / 1000;  // convert nanoseconds to microseconds
	totalTime += thistime;
	if (thistime > maxtime)
	    maxtime = thistime;
	if (thistime < mintime)
	    mintime = thistime;
      }

      System.out.println ("\nCPU time for " + iterations + " iterations: " + totalTime + " microseconds");
      System.out.println ("Average = " + (totalTime / iterations) + " microseconds per iteration");
      System.out.println ("Max = " + maxtime + " microseconds");
      System.out.println ("Min = " + mintime + " microseconds");
    }
  }
  
  public static void printByteArray (byte[] in) {
    System.out.print ("0x");
    for (int i = 0; i < in.length; i = i+2)
      System.out.printf ("%02x%02x ", in[i], in[i+1]);
  }


  public static long getCpuTime () {
  // Get CPU Time in nanoseconds

    ThreadMXBean bean = ManagementFactory.getThreadMXBean();
    return bean.isCurrentThreadCpuTimeSupported()?
      bean.getCurrentThreadCpuTime(): 0L;
  }


  public static long getUserTime () {
  // Get User Time in nanoseconds

    ThreadMXBean bean = ManagementFactory.getThreadMXBean();
    return bean.isCurrentThreadCpuTimeSupported()?
      bean.getCurrentThreadUserTime(): 0L;
  }


  public static long getSystemTime () {
  // Get System Time in nanoseconds

    ThreadMXBean bean = ManagementFactory.getThreadMXBean();
    return bean.isCurrentThreadCpuTimeSupported()?
      (bean.getCurrentThreadCpuTime() - bean.getCurrentThreadUserTime()): 0L;
  }
}
