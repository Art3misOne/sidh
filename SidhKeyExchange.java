
package sidh;

/**************************************************************************************************
 *
 * Implements supersingular isogeny Diffie-Hellman key exchange. Two constructors are provided. If
 * no arguments are passed, then it will default to using the parameters recommended by Microsoft
 * Research in the paper "Efficient algorithms for supersingular isogeny Diffie-Hellman" by 
 * Costello, Longa, and Naerhig. Alternate key exchange parameters can be used by passing them as
 * arguments to the constructor.
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.lang.System;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

class SidhKeyExchange {
  public static int PARTYA = 0;
  public static int PARTYB = 1;  

  int f;
  int lA;
  int lB;
  int eA;
  int eB;

  BigInteger prime;
  BigInteger orderB;
  BigInteger orderA;

  int obitsA;
  int obitsB;
    
  F2elm aGenPx;
  F2elm aGenQx;
  F2elm aGenDx;
  F2elm bGenPx;
  F2elm bGenQx;
  F2elm bGenDx;
    
  int maxA;
  int maxB;

  int maxIntPointsA;
  int maxIntPointsB;

  int[] splitsA;
  int[] splitsB;

  MontCurve baseCurve;


 public SidhKeyExchange() {
    setP503();
  }

    
  public SidhKeyExchange(String parameterID) {
    if (parameterID.equals ("sidhP751"))
      setP751 ();
    else if (parameterID.equals ("sidhP503"))
      setP503 ();
    else {
      System.out.println ("Unimplemented parameter set\n\n");
      System.exit(0);
    }
  }


  public void setP503() {
    // P503 parameters from SIKE NIST submission

    BigInteger x0, x1;

    f = 1;
    lA = 2;
    lB = 3;
    eA = 250;
    eB = 159;

    // prime = f*(lA^eA)*(lB^eB) - 1
    prime = new BigInteger("004066F541811E1E6045C6BDDA77A4D01B9BF6C87B7E7DAF13085BDA2211E7A0A" +
			   "BFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);

    try {
      Felm.setPrime (prime);
    } catch (InvalidFieldException ex) {
      System.out.println ("\nUnrecoverable error: Change in default parameters resulted in an " +
                            "invalid field.\n");
      System.exit (-99); 
    }

    orderA = new BigInteger ("400000000000000000000000000000000000000000000000000000000000000", 16);
    orderB = new BigInteger ("1019BD50604787981171AF769DE93406E6FDB21EDF9F6BC4C216F6888479E82B", 16);
    maxA = 125;
    maxB = 159;
    
    obitsA = 250;
    obitsB = 253;

    x0 = new BigInteger ("1f6d52a7563bb9356b98a116a0ca9775dbb7382eb29e24e45299d8939959eaeeb47" +
			 "ff3113f60882d12103e4b8b8cd2b97da14657ae8c128be82209d2ddfca9", 16);
    x1 = new BigInteger ("2d44c3fad24e4cbddc8a2d9de336a92a9912ee6d09e2dd5c33ab26d60a268ac91f3" +
			 "8e1af4c2d5bfa2b87dd55c8ca6019c6b0c08ed92b5aeb6c65a8e06e53e9", 16);
    aGenPx = new F2elm (x0, x1);

    x0 = new BigInteger ("97453912e12f3daf32eeffd618bd93d3bbbf399137bd39858cadefae382e42d6e60" +
			 "a62fd62417ad61a14b60db26125273ec980981325d86e55c45e3bb46b1", 16);
    aGenQx = new F2elm (x0, BigInteger.ZERO);

    x0 = new BigInteger ("173775ecbec79c78fd1ed5fe36075aace1f53f8ffb97d2a7e80dfc2875e77ec72d1" +
			 "d4a99e13353ec9d147badd96126948a72b30bdd7cebad7b54f8ddb5cd06", 16);
    x1 = new BigInteger ("2eaa224ddda149bbbb9089d2b2c471d068eca203465ce97dbc1c8ed0ebb0ff90e4f" +
			 "be7e266bba99cbae051797b4d35d28e36c1b1cb994aeeed1cb59fe5015", 16);
    aGenDx = new F2elm (x0, x1);

    x0 = new BigInteger ("21b7098b640a01d88708b729837e870cff9df6d4df86d86a7409f41156cb5f7b851" +
			 "4822730940c9b51e0d9821b0a67dd7ed98b9793685fa2e22d6d89d66a4e", 16);
    x1 = new BigInteger ("2f37f575bebbc33851f75b7ab5d89fc3f07e4df3cc52349804b8d17a17000a42fc6" +
			 "c5734b9fcfde669730f3e8569ceb53821d3e8012f7f391f57364f402909", 16);
    bGenPx = new F2elm (x0, x1);

    x0 = new BigInteger ("1e7d6ebceec9cfc47779affd696a88a971cdf3ec61e009df55caf4b6e01903b2cd1" +
			 "a12089c2ece106bdf745894c14d7e39b6997f70023e0a23b4b3787ef08f", 16);
    bGenQx = new F2elm (x0, BigInteger.ZERO);

    x0 = new BigInteger ("d4818d120a24abf48db51d129e6b1f24f4bbb2c16facc0c8c06323eeec2fa5b5e88" +
			 "7e17226417b1907310bfe6784fdebbac8c2a9abbe753f52259a7b7d70e", 16);
    x1 = new BigInteger ("19e75f0f03312d22cbbf153747525d89e5155babb8bf0c130cb567ca532f69aaf57" +
			 "ea7682b9957021d90414433abbeedc233e9082185781c16724c8c356777", 16);
    bGenDx = new F2elm (x0, x1);

    maxIntPointsA = 7;
    maxIntPointsB = 8;

    splitsA = new int[] {
      61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 
      4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 
      1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 
      1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 
      1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1 };

    splitsB = new int[] {
      71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 
      1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 
      5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 
      1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 
      2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 
      1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

    baseCurve = new MontCurve();
  }
    

  public void setP751() {
    // P751 parameters from SIKE NIST submission

    BigInteger x0, x1;

    f = 1;
    lA = 2;
    lB = 3;
    eA = 372;
    eB = 239;

    maxA = 186;
    maxB = 239;
    
    obitsA = 372;
    obitsB = 379;

    // prime = f*(lA^eA)*(lB^eB) - 1
    prime = new BigInteger("6fe5d541f71c0e12909f97badc668562b5045cb25748084e9867d6ebe876da959" +
                           "b1a13f7cc76e3ec968549f878a8eeafffffffffffffffffffffffffffffffffff" +
                           "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);

    try {
      Felm.setPrime (prime);
    } catch (InvalidFieldException ex) {
      System.out.println ("\nUnrecoverable error: Change in default parameters resulted in an " +
                            "invalid field.\n");
      System.exit (-99); 
    }

    // orderA = lA^eA, orderB = lB^eB
    orderA = new BigInteger ("1000000000000000000000000000000000000000000000000000000000000" +
                             "000000000000000000000000000000000", 16);
    orderB = new BigInteger ("6fe5d541f71c0e12909f97badc668562b5045cb25748084e9867d6ebe876d" +
                             "a959b1a13f7cc76e3ec968549f878a8eeb", 16);

    x0 = new BigInteger ("54921c31f0dc9531cb890fc5ec66df2e7f0d55761363c6e375da69b0682cabe5c" +
			 "0fffcbe6e1ad46563f042fa06b9f207fcf3cdd2673652828ff50c3f7b755c0be0" +
			 "72950d16ca747c146775c0267a401ffc738b03a49e9a36b39572afb363", 16);
    x1 = new BigInteger ("28849bc0d81e01993137a5b63d6e633c4e97ab4ff118ccf63dfe623092ac86b6d" +
			 "4a9b751797cba1a177500e9eb5af7852b7df02c334844d652efc4729178a1dbad" +
			 "8ca47bb7e757c6d43b799811a63bebe649c18101f03ad752cdcd73bf66", 16);
    
    aGenPx = new F2elm (x0, x1); 
    
    x0 = new BigInteger ("3e82027a38e9429c8d36ff46bcc93fa23f89f6be06d2b1317ad90438621783fdb" +
			 "7a4ad3e83e86cae096d5db822c98e561e008fa0e3f3b9ac2f40c56d6fa4a58a20" +
			 "449af1f1335661d14ab7347693632646086ce3acd54b0346f5cce233e9", 16);
    aGenQx = new F2elm (x0, BigInteger.ZERO);

    x0 = new BigInteger ("22a0b5a35a2b0c56135a7cec5cfb97964a7c6226fe909f374362a8eca3ab14a1b" +
			 "7b0c87ac875dce5888d83b623bf0011a4ac138f62ef6b2d2d84f636548a9f920f" +
			 "238336e5a36e45e4055940e3c94385b8fc5374396432eef2ae178cefdd", 16);
    x1 = new BigInteger ("f9c4afcda809c3358b096b250c69b20310fdf2ef631711aa4efec49a4e76483f3" +
			 "20b793f2ebc63365eed14aa3f6ea33feb56796f011ba6c6dfb4d0a00aac4d2786" +
			 "646d914ad026cbb4a592ec74b5485372e51382d44528dd491b83d9547", 16);
    aGenDx = new F2elm (x0, x1);
    
    x0 = new BigInteger ("5fd1a3c4dd0f630974196fed3519152bc7098b9e2b121eca46bd10a5cc9f4bcc6" +
			 "c689b8e4c063b3798075fcee6edaa9eb108b3cd00495cf04dd8ce4a08fbe685a1" +
			 "27d40e45f4cf45098a578deb44368699394c43bfc9bc5e00052f78e8d", 16);
    x1 = new BigInteger ("2b88a03360b3389547732c9140c05dea6516881fe108211be887cc43fcb80c06a" +
			 "1d86ff5457d3bb7db936394ec33821aa39333a60af84b537974cfa0ba8287d699" +
			 "d2bf79ba559026c64a6ed610501d2357c10b9a6c8f837424922275acbf", 16);
    bGenPx = new F2elm (x0, x1);

    x0 = new BigInteger ("2f1d80ef06ef960a01ab8ff409a2f8d5bce859ed725de145fe2d525160e0a3ad8" +
			 "e17b9f9238cd5e69cf26df237429bd3778659023b9ecb610e30288a7770d3785a" +
			 "aaa4d646c576aecb94b919aeedd9e1df566c1d26d376ed2325dcc93103", 16);
    bGenQx = new F2elm (x0, BigInteger.ZERO);

    x0 = new BigInteger ("77b3bb69009428a327d43ca60169715f547454f88cd017b32df58a7252c2b3c3d" +
			 "00d52ccd3133d54041d8bcaea291f2057202328712cd395575cd7ccd3ce70c0a1" +
			 "ebf633ba946559458878f41f9fdd1727e2c31125b2fe5b71306704829", 16);
    x1 = new BigInteger ("6d91393a57dbf47fd6dcf841f17ecd719cae1d33c6832a75b0f168855bcc38d2a" +
			 "4792dff9bc86deaca10b1aa808d539b167d73bba32168687fa3f85ae93a1adde5" +
			 "bd1fd5b681dcc6c34454d4496976c22d80c95e42b12576fc0fb4074b9f", 16);
    bGenDx = new F2elm (x0, x1);
    
    maxIntPointsA = 8;
    maxIntPointsB = 10;

    splitsA = new int[] { 
	80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
	1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
	1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 
	1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 
	33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 
	1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 
	1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

    splitsB = new int[] { 
	112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
	1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
	1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
	1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 
	2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 
	1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 
	15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 
	1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
	2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

    baseCurve = new MontCurve();
  }


  public MontCurve getCurve() {
    return baseCurve;
  }


  public BigInteger getOrderA () {
    return orderA;
  }


  public BigInteger getOrderB () {
    return orderB;
  }

    
  public F2elm[] getGenA() {
    F2elm genA[] = new F2elm[] { aGenPx, aGenQx, aGenDx };
    return genA;
  }


  public F2elm[] getGenB() {
    F2elm genB[] = new F2elm[] { bGenPx, bGenQx, bGenDx };
    return genB;
  }
    

  public int getMIPA() {
    return maxIntPointsA;
  }


  public int getMIPB() {
    return maxIntPointsB;
  }


  public int getMaxA() {
    return maxA;
  }


  public int getMaxB() {
    return maxB;
  }


  public int[] getSplitsA() {
    return splitsA;
  }


  public int[] getSplitsB() {
    return splitsB;
  }


  public int getObitsA() {
    return obitsA;
  }


  public int getObitsB() {
    return obitsB;
  }


  public SidhKeyPair generateKeyPair (int aOrB) {
    return new SidhKeyPair (aOrB, this);
  }


  public SidhKeyPair generateKeyPair (int aOrB, SidhPrivateKey prKey) {
    return new SidhKeyPair (aOrB, prKey, this);
  }


  public byte[] calculateAgreementA (SidhPrivateKey privKeyA, SidhPublicKey pubKeyB) {
    int i, ii = 0, row, m, index = 0, ptsIdx[], npts = 0;
    F2Point r, pts[];
    F2elm aB, pkB0, pkB1, pkB2, two;
    FourIsogeny fourIsog;

    pkB0 = new F2elm (pubKeyB.getP ());
    pkB1 = new F2elm (pubKeyB.getQ ());
    pkB2 = new F2elm (pubKeyB.getD ());

    aB = MontCurve.recoverA (pkB0, pkB1, pkB2);

    fourIsog = new FourIsogeny (aB);
    fourIsog.updateA24 ();
    two = F2elm.leftShift (F2elm.ONE, 1);              
    fourIsog.setA24plus (F2elm.add (two, aB));         
    fourIsog.setC24 (F2elm.leftShift (two, 1));        
    
    pts = new F2Point[maxIntPointsA];
    ptsIdx = new int[maxIntPointsA];

    r = fourIsog.ladder3pt (pkB0, pkB1, pkB2, privKeyA.getKey (), obitsA);

    for (row = 1; row < maxA; row++) {
      while (index < maxA - row) {
        pts[npts] = r;
        ptsIdx[npts++] = index;
        m = splitsA[ii++];
        r = fourIsog.xDble (r, 2*m);
        index += m;
      }

      fourIsog.get4Isog (r);
      
      for (i = 0; i < npts; i++) 
        pts[i] = fourIsog.eval4Isog (pts[i]);

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    fourIsog.get4Isog (r);
    fourIsog.updateAC (4);

    return fourIsog.jInv().toByteArray();
  }


  public byte[] calculateAgreementB (SidhPrivateKey privKeyB, SidhPublicKey pubKeyA) {
    int i, ii = 0, row, m, index = 0, ptsIdx[], npts = 0;
    F2Point r, pts[];
    F2elm pkA0, pkA1, pkA2, aA, temp;
    ThreeIsogeny threeIsog;

    pts = new F2Point[maxIntPointsB];
    ptsIdx = new int[maxIntPointsB];

    pkA0 = pubKeyA.getP ();
    pkA1 = pubKeyA.getQ ();
    pkA2 = pubKeyA.getD ();

    aA = MontCurve.recoverA (pkA0, pkA1, pkA2);
    
    threeIsog = new ThreeIsogeny (aA);
    threeIsog.updateA24 ();
    temp = F2elm.leftShift (F2elm.ONE, 1);
    threeIsog.setA24plus (F2elm.add (aA, temp));
    threeIsog.setA24minus (F2elm.sub (aA, temp));
    
    r = threeIsog.ladder3pt(pkA0, pkA1, pkA2, privKeyB.getKey (), obitsB);

    for (row = 1; row < maxB; row++) {
      while (index < maxB - row) {
        pts[npts] = r;
        ptsIdx[npts++] = index;
        m = splitsB[ii++];
        r = threeIsog.xTple (r, m);
        index += m;
      }

      threeIsog.get3Isog (r);

      for (i = 0; i < npts; i++) 
        pts[i] = threeIsog.eval3Isog (pts[i]);

      r = pts[npts-1];
      index = ptsIdx[npts-1];
      npts--;
    }

    threeIsog.get3Isog (r);
    threeIsog.updateAC (3);

    return threeIsog.jInv().toByteArray();
  }

    
  public static void writeKeyToFile (String filename, SidhPublicKey PubKey) {
    OutputStream fStream = null;

    try {
      byte[] keyBytes = PubKey.serialize();
      File fh = new File (filename);
      if (!fh.exists()) {
	fh.createNewFile ();
      }
      
      fStream = new FileOutputStream (fh);
      fStream.write (keyBytes);
      fStream.flush();
    } catch (IOException e) {
      e.printStackTrace();
    } finally {
      try {
	if (fStream != null) fStream.close();
      } catch (IOException e) {
      }
    }
  }
}
