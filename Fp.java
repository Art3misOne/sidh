
package sidh;

/**************************************************************************************************
 *
 * Implements elements over finite field GF(p) and the quadratic extension GF(p^2) 
 *  
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;
import java.lang.System;
import java.lang.Math;
  import java.security.SecureRandom;

class InvalidFieldException extends Exception {
}


class Felm {
  /* Elements of GF(p) */

  // Set p to a placeholder value until the prime has been set
  public static BigInteger p = BigInteger.valueOf(2);
  private BigInteger value;
  public static int primesize = 1;

  public static final Felm ZERO = new Felm (BigInteger.ZERO);
  public static final Felm ONE = new Felm (BigInteger.ONE);
  
    
  public Felm (BigInteger v) {
    value = v.mod(p);
  }


  public Felm (Felm a) {
    value = a.fpGetValue();
  }


  public Felm (long v) {
    value = BigInteger.valueOf(v).mod(p);
  }
    

  public Felm (byte[] bytes) {
    value = new BigInteger (bytes);
  }


  public Felm (SecureRandom rnd) {
    value = genRandom (p);
  }
    

  public static void setPrime (BigInteger pr) throws InvalidFieldException {
    if (pr.isProbablePrime(10))
      p = pr;
    else
      throw new InvalidFieldException();

    primesize = (pr.bitLength() / 8) + 1; 
  }


  public static BigInteger getPrime() {
    return p;
  }


  public BigInteger fpGetValue() {
    return value;
  }


  public static BigInteger genRandom (BigInteger bound) {
    // Gen random values up to the same bit length as the bound until the value generated is
    // strictly less than the bound. Average expected number of calls is less than 2.
    
    SecureRandom rnd = new SecureRandom();
    int numBits = bound.bitLength();
    
    BigInteger randval = new BigInteger (numBits, rnd);
        
    while (randval.compareTo(bound) >= 0)
        randval = new BigInteger (numBits, rnd);

    return randval;
  }


  public void randomize () {
    value = genRandom (p);
  }
    
    
  public static Felm add (Felm x, Felm y) {
    Felm z = new Felm (x);
    z.fpAddInPlace (y);
    return z;
  }


  public void fpAddInPlace (Felm y) {
    value = value.add (y.value);
    value = value.mod (p);
  }


  public static Felm sub (Felm x, Felm y) {
    Felm z = new Felm (x);
    z.fpSubInPlace (y);
    return z;
  }

    
  public void fpSubInPlace (Felm y) {
    value = value.add (p).subtract (y.value);
    value = value.mod (p);
  }


  public static Felm mult (Felm x, Felm y) {
    Felm z = new Felm (x);
    z.fpMultInPlace (y);
    return z;
  }


  public void fpMultInPlace (Felm y) {
    value = value.multiply (y.value);
    value = value.mod (p);
  }

    
  public static Felm sqr (Felm x) {
    return mult (x, x);
  }


  public void fpSqrInPlace () {
    fpMultInPlace (this);
  }
    

  public boolean fpIsZero() {
    return value.equals (BigInteger.ZERO);
  }


  public boolean fpIsEven() {
    BigInteger c = value.and (BigInteger.ONE);
    return c.equals (BigInteger.ZERO);
  }


  public boolean fpIsOdd() {
    return !fpIsEven();
  }


  public boolean fpEquals (Felm y) {
    return value.equals (y.value);
  }


  public boolean fpIsLessThan (Felm y) {
    return value.compareTo (y.value) == -1;
  }


  public boolean fpIsGreaterThan (Felm y) {
    return value.compareTo (y.value) == 1;
  }


  public static Felm negate (Felm x) {
    return new Felm (p.subtract (x.value));
  }


  public void fpNegateInPlace () {
    value = p.subtract (value);
  }
    

  public static Felm inverse (Felm x) {
    Felm z = new Felm (x);
    z.fpInverseInPlace ();
    return z;
  }


  public void fpInverseInPlace () {
    value = value.modInverse (p);
  }
    

  public static Felm div2 (Felm x) {
    Felm z = new Felm (x);
    z.fpDiv2InPlace ();
    return z;
  }


  public void fpDiv2InPlace () {
    if (fpIsOdd ())
      value = value.add (p);
    value = value.shiftRight (1);    
  }


  public static Felm leftShift (Felm x, int shiftBy) {
    Felm z = new Felm (x);
    z.fpLeftShiftInPlace (shiftBy);
    return z;
  }


  public void fpLeftShiftInPlace (int shiftBy) {
    value = value.shiftLeft (shiftBy);
  }
    

  public static Felm rightShift (Felm x, int shiftBy) {
    Felm z = new Felm (x);
    z.fpRightShiftInPlace (shiftBy);
    return z;
  }


  public void fpRightShiftInPlace (int shiftBy) {
    value = value.shiftRight (shiftBy);
  }
    

  public Felm fpSwap (Felm y, BigInteger option) {
    // Constant time swap regardless of whether option is 0 or 1
    BigInteger temp, mask, yval;

    yval = y.value;
    mask = option.negate();                  // option = 1 => mask = 1...1 because BigInteger
                                             // automatically sign extends as necessary

    temp = mask.and(value.xor(yval));        // temp = mask & (this.value xor y.value)
    value = temp.xor(value);
    yval = temp.xor(yval);

    return new Felm (yval);
  } 
    

  public String toString() {
    return "0x" + value.toString(16);
  }


  public byte[] toByteArray() {
    // Returns the same size array regardless of the value. Zero pad the highbits.
    byte[] retval = new byte[primesize];
    Arrays.fill (retval, (byte) 0);                    

    int eltsize = (value.bitLength() / 8) + 1; 
    int offset = primesize - eltsize;

    System.arraycopy (value.toByteArray(), 0, retval, offset, eltsize);

    return retval;
  }
} 


class F2elm {
  /* Elements of the quadratic extension field GF(p^2): x0 + x1*i */

  private Felm x0;
  private Felm x1;

  public static final F2elm ZERO = new F2elm (Felm.ZERO, Felm.ZERO);
  public static final F2elm ONE = new F2elm (Felm.ONE, Felm.ZERO);


  public F2elm (BigInteger a0, BigInteger a1) {
    x0 = new Felm (a0);
    x1 = new Felm (a1);
  }


  public F2elm (Felm a0, Felm a1) {
    x0 = new Felm (a0);
    x1 = new Felm (a1);
  }


  public F2elm (F2elm a) {
    x0 = new Felm (a.x0);
    x1 = new Felm (a.x1);
  }


  public F2elm (long v0, long v1) {
    x0 = new Felm (v0);
    x1 = new Felm (v1);
  }
    

  public F2elm (byte[] bytes) {
    int len = (bytes.length) / 2;
    x0 = new Felm (Arrays.copyOfRange (bytes, 0, len));
    x1 = new Felm (Arrays.copyOfRange (bytes, len, 2*len));
  }


  public F2elm (SecureRandom rnd) {
    x0 = new Felm (rnd);
    x1 = new Felm (rnd);
  }
    

  public Felm f2Get0() {
    return x0;
  }


  public Felm f2Get1() {
    return x1;
  }


  public boolean f2Equals (F2elm y) {
    return x0.fpEquals(y.x0) && x1.fpEquals(y.x1);
  }


  public static F2elm add (F2elm x, F2elm y) {
    F2elm z = new F2elm (x);
    z.f2AddInPlace (y);
    return z;
  }


  public void f2AddInPlace (F2elm y) {
    x0.fpAddInPlace (y.x0);
    x1.fpAddInPlace (y.x1);
  }

    
  public static F2elm sub (F2elm x, F2elm y) {
    F2elm z = new F2elm (x);
    z.f2SubInPlace (y);
    return z;
  }


  public void f2SubInPlace (F2elm y) {
    x0.fpSubInPlace (y.x0);
    x1.fpSubInPlace (y.x1);
  }

    
  public static F2elm negate (F2elm x) {
    F2elm y = new F2elm (x);
    y.f2NegateInPlace ();
    return y;
  }


  public void f2NegateInPlace () {
    x0.fpNegateInPlace ();
    x1.fpNegateInPlace ();
  }
    

  public static F2elm sqr (F2elm x) {
    F2elm y = new F2elm (x);
    y.f2SqrInPlace ();
    return y;
  }


  public void f2SqrInPlace () {
    Felm t1, t2, t3;

    t1 = Felm.add (x0, x1);                       // t1 = x0 + x1
    t2 = Felm.sub (x0, x1);                       // t2 = x0 - x1
    t3 = Felm.leftShift (x0, 1);                  // t3 = 2 * x0

    x0 = Felm.mult (t1, t2);                      // x0 = (x0+x1)(x0-x1)
    x1.fpMultInPlace (t3);                        // x1 = 2*x0*x1
  }
    

  public static F2elm mult (F2elm y, F2elm z) {
    F2elm x = new F2elm (y);
    x.f2MultInPlace (z);
    return x;
  }


  public void f2MultInPlace (F2elm y) {
    // compute c = this * y = (x0 + i*x1) * (y0 + i*y1) = x0y0 - x1y1 + i*(x0y1 + x1y0)
    Felm t1, t2, c0, c1, y0, y1;

    y0 = y.x0;
    y1 = y.x1;
    
    t1 = Felm.mult (x0, y0);                      
    t2 = Felm.mult (x1, y1);                      
    c0 = Felm.sub (t1, t2);                       

    // Using extra additions, but fewer multiplications
    t1.fpAddInPlace (t2);
    t2 = Felm.add (x0, x1);
    x1 = Felm.add (y0, y1);
    x1.fpMultInPlace (t2);
    x1.fpSubInPlace (t1);
    x0 = c0;
  }
    

  public static F2elm rightShift (F2elm y, int n) {
    F2elm x = new F2elm (y);
    x.f2RightShiftInPlace (n);
    return x;
  }


  public void f2RightShiftInPlace (int n) {
    x0.fpRightShiftInPlace (n);
    x1.fpRightShiftInPlace (n);
  }
    

  public static F2elm leftShift (F2elm y, int n) {
    F2elm x = new F2elm (y);
    x.f2LeftShiftInPlace (n);
    return x;
  }


  public void f2LeftShiftInPlace (int n) {
    x0.fpLeftShiftInPlace (n);
    x1.fpLeftShiftInPlace (n);
  } 
    

  public boolean f2IsEven () {
    return x0.fpIsEven() && x1.fpIsEven();
  }


  public static F2elm div2 (F2elm y) {
    F2elm x = new F2elm (y);
    x.f2Div2InPlace ();
    return x;
  }


  public void f2Div2InPlace () {
    x0.fpDiv2InPlace ();
    x1.fpDiv2InPlace ();
  }
    

  public static F2elm inverse (F2elm y) {
    F2elm x = new F2elm (y);
    x.f2InverseInPlace ();
    return x;
  }


  public void f2InverseInPlace () {
    Felm t0, t1;

    t0 = Felm.sqr (x0);
    t1 = Felm.sqr (x1);
    t0.fpAddInPlace (t1);
    t0.fpInverseInPlace ();
    x1.fpNegateInPlace ();
    x0.fpMultInPlace (t0);
    x1.fpMultInPlace (t0);
  }
    

  public static F2elm[] inv3Way (F2elm z0, F2elm z1, F2elm z2) {
    // Compute simultaneous inversion of 3 elements

    F2elm t0, res[] = new F2elm[3];

    t0 = mult (z0, z1);                    // t0 = z0*z1 
    res[1] = mult (t0, z2);                // res1 = z0*z1*z2
    res[2] = inverse (res[1]);             // res2 = 1/(z0*z1*z2)
    res[1] = mult (res[2], z2);            // res1 = 1/(z0*z1)
    res[0] = mult (res[1], z1);            // res0 = 1/z0
    res[1].f2MultInPlace (z0);             // res1 = 1/z1
    res[2].f2MultInPlace (t0);             // res2 = 1/z2

    return res;
  }
    

  public static F2elm[] inv4Way (F2elm z0, F2elm z1, F2elm z2, F2elm z3) {
    // Compute simultaneous inversion of 4 elements
    
    F2elm res[] = new F2elm[4];
    
    res[0] = mult(z0, z1);                 // res0 = z0*z1
    res[1] = mult(z2, z3);                 // res1 = z2*z3
    res[2] = mult(res[0], res[1]);         // res2 = z0*z1*z2*z3
    res[3] = inverse(res[2]);              // res3 = 1/(z0*z1*z2*z3)
    res[2] = mult(res[1], res[3]);         // res2 = 1/(z2*z3)
    res[3].f2MultInPlace(res[0]);          // res3 = 1/(z0*z1)
    res[0] = mult(res[2], z1);             // res0 = 1/z0
    res[1] = mult(res[2], z0);             // res1 = 1/z1
    res[2] = mult(res[3], z3);             // res2 = 1/z2
    res[3].f2MultInPlace(z2);              // res3 = 1/z3

    return res;
  }


  public F2elm f2Swap (F2elm y, BigInteger option) {
    // Constant time swap regardless of whether option is 0 or 1
    Felm y0, y1;

    y0 = x0.fpSwap(y.x0, option);
    y1 = x1.fpSwap(y.x1, option);

    return new F2elm (y0, y1);
  }    


  public static F2elm select (F2elm x, F2elm y, BigInteger option) {
    // Return x if option = 0 and y if option = 1
    BigInteger y0, y1, z0, z1, mask, bix0, bix1;

    mask = option.negate ();                // if option = 1 then mask = 1...1 because BigInteger
                                            // automatically sign extends as necessary

    // Get x0, x1, y0, y1 as their BigInteger representations
    y0 = (y.x0).fpGetValue ();
    y1 = (y.x1).fpGetValue ();
    bix0 = (x.x0).fpGetValue ();                 
    bix1 = (x.x1).fpGetValue ();                 

    z0 = bix0.xor (y0);                     // z0 = x0 xor y0
    z1 = bix1.xor (y1);                     // z1 = x1 xor y1

    z0 = z0.and (mask);                     // if mask = 0 then z0 = 0 else z0 = x0 xor y0 
    z1 = z1.and (mask);                     // if mask = 0 then z1 = 0 else z1 = x1 xor y1 

    z0 = bix0.xor (z0);                     // if mask = 0 then z0 = x0 
                                            // else z0 = x0 xor x0 xor y0 = y0

    z1 = bix1.xor (z1);                     // if mask = 0 then z1 = x1 
                                            // else z1 = x1 xor x1 xor y1 = y1

    return new F2elm (z0, z1);
  }


  public String toString() {
    return "[" + x0 + "," + x1 + "]";
    //return x1 + "*i + " + x0;
  }


  public byte[] toByteArray() {
    byte[] retval = new byte[2*Felm.primesize];
    System.arraycopy (x0.toByteArray(), 0, retval, 0, Felm.primesize);
    System.arraycopy (x1.toByteArray(), 0, retval, Felm.primesize, Felm.primesize);
    return retval;
  }
}


