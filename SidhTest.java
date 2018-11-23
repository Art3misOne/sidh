
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
  public static boolean debug = false;
    
  public static void main (String[] args) {
    // Using default parameters
    SidhKeyExchange kex = new SidhKeyExchange("sidhP503");
    SidhKeyPair keysA, keysB;
    SidhPublicKey reconstructedApub, publicB;
    SidhPrivateKey reconstructedApriv, privateB;
    byte[] sharedA, sharedB, pubKeyBytes, privKeyBytes, knownAnswer;
    
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
        System.out.println ("\nPublic key reconstruction successful");
      else
        System.out.println ("\nPublic key reconstruction unsuccessful");

      if (reconstructedApriv.privateKeyEquals (keysA.getPrivateKey()))
        System.out.println ("Private key reconstruction successful\n");
      else
        System.out.println ("Private key reconstruction unsuccessful\n");      
      
      sharedA = kex.calculateAgreementA (keysA.getPrivateKey(), keysB.getPublicKey());
      sharedB = kex.calculateAgreementB (keysB.getPrivateKey(), keysA.getPublicKey());

      knownAnswer = new byte[] {
	(byte) 0x27, (byte) 0x98, (byte) 0x02, (byte) 0xa3, (byte) 0xda, (byte) 0xd9, (byte) 0x05,
	(byte) 0xf7, (byte) 0xc7, (byte) 0x03, (byte) 0x88, (byte) 0xa9, (byte) 0x38, (byte) 0x92,
	(byte) 0x0a, (byte) 0x77, (byte) 0xc0, (byte) 0x7e, (byte) 0x23, (byte) 0x21, (byte) 0x36,
	(byte) 0x25, (byte) 0xc3, (byte) 0xd7, (byte) 0x3d, (byte) 0x5f, (byte) 0x92, (byte) 0x11,
	(byte) 0x6c, (byte) 0x3f, (byte) 0x3f, (byte) 0x75, (byte) 0xe4, (byte) 0x21, (byte) 0x3d,
	(byte) 0x69, (byte) 0xa7, (byte) 0x85, (byte) 0xe6, (byte) 0x7b, (byte) 0xa0, (byte) 0x02,
	(byte) 0xce, (byte) 0x03, (byte) 0xca, (byte) 0x83, (byte) 0xbb, (byte) 0x25, (byte) 0x6c,
	(byte) 0xa1, (byte) 0xc2, (byte) 0xaa, (byte) 0x65, (byte) 0x15, (byte) 0x8b, (byte) 0x9c,
	(byte) 0x2f, (byte) 0x33, (byte) 0xff, (byte) 0x39, (byte) 0x78, (byte) 0xe2, (byte) 0x50,
	(byte) 0x25, (byte) 0xb5, (byte) 0xa9, (byte) 0x7d, (byte) 0x2a, (byte) 0x07, (byte) 0xba,
	(byte) 0x96, (byte) 0xb3, (byte) 0x1e, (byte) 0x58, (byte) 0xb7, (byte) 0xd5, (byte) 0x40,
	(byte) 0xc5, (byte) 0x4a, (byte) 0x28, (byte) 0xd5, (byte) 0xba, (byte) 0x6e, (byte) 0xd2,
	(byte) 0x3b, (byte) 0x86, (byte) 0x0f, (byte) 0x04, (byte) 0x2a, (byte) 0x35, (byte) 0x2e,
	(byte) 0x79, (byte) 0xa7, (byte) 0x92, (byte) 0x90, (byte) 0xd1, (byte) 0xf9, (byte) 0xe2,
	(byte) 0xcf, (byte) 0xb6, (byte) 0xa1, (byte) 0xd4, (byte) 0x31, (byte) 0x5a, (byte) 0xf0,
	(byte) 0xe9, (byte) 0x4d, (byte) 0x48, (byte) 0x57, (byte) 0xf8, (byte) 0x71, (byte) 0xbd,
	(byte) 0xa7, (byte) 0xf6, (byte) 0xd2, (byte) 0xc2, (byte) 0xdb, (byte) 0xb7, (byte) 0x32,
	(byte) 0xbe, (byte) 0x6b, (byte) 0x01, (byte) 0x88, (byte) 0x98, (byte) 0x4e, (byte) 0x47
      };
      
      if (Arrays.equals (sharedA, sharedB)) {
	if (Arrays.equals (sharedA, knownAnswer))
	    System.out.println ("\nShared secrets match between users and match known answer :)\n");
	else {
	  System.out.print ("\nShared secrets match, but not against known answer =/\n");
	  if (debug) {
	    System.out.print ("Known answer = ");
	    printByteArray (knownAnswer);
	    System.out.println ("\n");
	  }
	}
	
        if (debug) {
	  System.out.print ("Shared secret = ");
	  printByteArray (sharedA);
	  System.out.println ("\n");
	}
      }

      else {
	System.out.println ("Shared secrets do not match :(\n");
	if (debug) {
	  System.out.print ("SharedA = ");
	  printByteArray (sharedA);
	  System.out.print ("\nSharedB = ");
	  printByteArray (sharedB);
	  System.out.println ("\n");
	}
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
      System.out.println ("Min = " + mintime + " microseconds\n");
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
