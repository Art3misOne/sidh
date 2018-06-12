
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
  public static int primesize = 1;           // bytes needed to represent the prime

  // Frequently used constants
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

    
  public Felm fpAdd (Felm y) {
    return new Felm (value.add (y.value));
  }


  public Felm fpSub (Felm y) {
    return fpAdd (y.fpNegate());             
  }


  public Felm fpMult (Felm y) {
    return new Felm (value.multiply (y.value));
  }


  public Felm fpSqr() {
    return fpMult(this);
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


  public Felm fpNegate() {
    return new Felm (p.subtract(value));
  }


  public Felm fpInverse() {
    return new Felm (value.modInverse(p));
  }


  public Felm fpDiv2 () {
    BigInteger v = value;

    if (fpIsOdd()) 
      v = v.add(p);                          // If this is odd, adding p makes it even
    v = v.shiftRight(1);                     // Divide by two = right shift by one

    return new Felm (v);
  }


  public Felm fpLeftShift (int shiftBy) {
    BigInteger v = value;
    v = v.shiftLeft (shiftBy);
    return new Felm (v);
  }


  public Felm fpRightShift (int shiftBy) {
    BigInteger v = value;
    v = v.shiftRight (shiftBy);
    return new Felm (v);
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


  public Felm fpSqrt () {
    // Compute square root using Tonelli-Shanks algorithm
    Felm c, z, t, r, b, tsq;
    BigInteger q, exp, exp2;
    int s, m, i;

    if (value.equals (BigInteger.ZERO))
      return Felm.ZERO;
    
    if (isQuadraticResidue() == false)
      return null;

    q = p.subtract (BigInteger.ONE);
    s = numPowersOf2 (q);                              // s = number of times 2 divides (p-1)
    q = q.shiftRight (s);                              // q = (p-1) / (2^s)

    z = quadraticNonResidue ();

    exp = q.add (BigInteger.ONE);
    exp = exp.shiftRight (1);                          // exp = (q+1)/2
    
    m = s;
    c = z.fpPow (q);                                   // c = z^q
    t = fpPow (q);                                     // t = value^q
    r = fpPow (exp);                                   // r = value^((q+1)/2)

    while (t.fpEquals (Felm.ONE) == false) {
      i = t.findLog2Order ();
      
      exp = BigInteger.ONE.shiftLeft (m-i-1);          // exp = 2^(m-i-1)
      b = c.fpPow (exp);                               // b = c^(2^(m-i-1))
      
      m = i;
      c = b.fpSqr ();                                  // c = b^2
      t = t.fpMult (c);                                // t = t*b^2
      r = r.fpMult (b);                                // r = r*b
    }

    return r;
  }


  private boolean isQuadraticResidue () {
    Felm c;
    BigInteger d;

    d = p.subtract (BigInteger.ONE);
    d = d.shiftRight (1);
    c = fpPow (d);

    if (c.fpEquals (Felm.ONE))
      return true;
    return false;
  }
    

  private int numPowersOf2 (BigInteger n) {
    BigInteger q = n;
    int s = 0;
    
    while (q.testBit(0) == false) {                    // while q is even
      s = s + 1;                                       
      q = q.shiftRight (1);                            // q = q/2
    }

    return s;
  }


  private Felm quadraticNonResidue () {
    Felm z;
    SecureRandom rnd = new SecureRandom ();

    do {
      z = new Felm (rnd);                              // z <- random
    } while (z.fpEquals (Felm.ZERO) || z.isQuadraticResidue ());

    return z;
  }


  private int findLog2Order () {
    Felm tsq;
    int i = 0;

    tsq = new Felm (this);
    while (tsq.fpEquals (Felm.ONE) == false) {
      i = i + 1;
      tsq = tsq.fpSqr ();
    }

    return i;
  }
    

  public Felm fpPow (BigInteger exp) {
    Felm result = Felm.ONE;
    Felm base = new Felm (this);
    BigInteger e = exp;

    while (e.compareTo (BigInteger.ZERO) == 1) {
      if (e.testBit (0) == true)
	result = result.fpMult (base);
      e = e.shiftRight (1);
      base = base.fpSqr ();
    }

    return result;
  }

    
  public String toString() {
    return value.toString();
  }


  // Creates a byte[] representation of the field element. Returns the same size array regardless
  // of the element so each one is the same number of bytes as needed to represent the prime.
  public byte[] toByteArray() {
    byte[] retval = new byte[primesize];
    Arrays.fill (retval, (byte) 0);                    

    int eltsize = (value.bitLength() / 8) + 1; 
    int offset = primesize - eltsize;

    // The return array was filled with zeros so copy from offset to the end of the element. This
    // way the zero padding is at the beginning and not the end since adding zeros to the end of a
    // binary representation of a number multiplies the number by a power of 2.
    System.arraycopy (value.toByteArray(), 0, retval, offset, eltsize);

    return retval;
  }
} 


class F2elm {
  /* Elements of the quadratic extension field GF(p^2): x0 + x1*i */

  private final Felm x0;
  private final Felm x1;

  // Constants of the class
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


  public F2elm f2Add (F2elm y) {
    Felm newx0, newx1;

    newx0 = x0.fpAdd(y.x0);
    newx1 = x1.fpAdd(y.x1);

    return new F2elm (newx0, newx1);
  }


  public F2elm f2Sub (F2elm y) {
    Felm newx0, newx1;

    newx0 = x0.fpSub(y.x0);
    newx1 = x1.fpSub(y.x1);

    return new F2elm (newx0, newx1);
  }


  public F2elm f2Negate() {
    return new F2elm (x0.fpNegate(), x1.fpNegate());
  }


  public F2elm f2Sqr () {
    // Compute c = this * this = (x0 + i*x1)^2 = x0^2 - x1^2 + i*2x0x1
    Felm t1, t2, t3, c0, c1;

    t1 = x0.fpAdd(x1);                       // t1 = x0 + x1
    t2 = x0.fpSub(x1);                       // t2 = x0 - x1
    t3 = x0.fpAdd(x0);                       // t3 = 2 * x0
    c0 = t1.fpMult(t2);                      // c0 = (x0+x1)(x0-x1)
    c1 = t3.fpMult(x1);                      // c1 = 2*x0*x1

    return new F2elm (c0, c1);
  }
    

  public F2elm f2Mult (F2elm y) {
    // compute c = this * y = (x0 + i*x1) * (y0 + i*y1) = x0y0 - x1y1 + i*(x0y1 + x1y0)
    Felm t1, t2, c0, c1, y0, y1;

    y0 = y.x0;
    y1 = y.x1;

    t1 = x0.fpMult(y0);                      // t1 = x0 * y0
    t2 = x1.fpMult(y1);                      // t2 = x1 * y1
    c0 = t1.fpSub(t2);                       // c0 = (x0*y0) - (x1*y1))

    // Using extra additions, but fewer multiplications
    t1 = t1.fpAdd(t2);                       // t1 = (x0*y0) + (x1*y1)
    t2 = x0.fpAdd(x1);                       // t2 = x0 + x1
    c1 = y0.fpAdd(y1);                       // c1 = y0 + y1
    c1 = t2.fpMult(c1);                      // c1 = (x0+x1)(y0+y1)
    c1 = c1.fpSub(t1);                       // c1 = (x0+x1)(y0+y1) - x0*y0 - x1*y1
                                             //      = (x0*y1) + (x1*y0)

    return new F2elm (c0, c1);
  }


  public F2elm f2MultByi () {
    return new F2elm (x1.fpNegate(), x0);
  }
    

  public boolean f2IsEven () {
    return x0.fpIsEven() && x1.fpIsEven();
  }


  public F2elm f2Div2 () {
    Felm c0, c1;

    c0 = x0.fpDiv2();
    c1 = x1.fpDiv2();

    return new F2elm (c0, c1);
  }


  public F2elm f2Sqrt () {
    BigInteger pmod4, three;
    
    if (this.f2Equals (F2elm.ZERO))
      return F2elm.ZERO;

    if (this.isQuadraticResidue () == false)
      return null;

    three = BigInteger.valueOf (3);
    pmod4 = Felm.getPrime ().and (three);
    
    if (pmod4.equals (three))
      return f2Sqrt3Mod4 ();
    else
      return f2Sqrt1Mod4 ();
  }


  private F2elm f2Sqrt3Mod4 () {
    F2elm a1, a2, alpha, b, neg1;
    BigInteger d, exp;

    neg1 = F2elm.ONE.f2Negate();

    d = Felm.getPrime().shiftRight (2);    // d = floor (p/4) = (p-3)/4
    a1 = f2Pow (d);                        // a1 = x^d = (x0,x1)^d
    a2 = f2Mult (a1);                      // a2 = x*a1 = x^(d+1) = x^((p+1)/4)
    alpha = a1.f2Mult (a2);                // alpha = a1*a2 = x^(2d+1) = x^((p-1)/2)

    if (alpha.f2Equals (neg1))
      return a2.f2MultByi ();              // return a2*i = x^(d+1) * i

    exp = Felm.getPrime ().shiftRight (1); // exp = (p-1)/2
    b = alpha.f2Add (F2elm.ONE);           // b = alpha + 1
    b = b.f2Pow (exp);                     // b = b^exp
    return b.f2Mult (a2);                  // return b * a2
  }


  private F2elm f2Sqrt1Mod4 () {
    F2elm a;
    Felm a0, a1;

    a0 = x0.fpSqr ();                      // a0 = x0^2 
    a1 = x1.fpSqr ();                      // a1 = x1^2
    a1 = a0.fpAdd (a1);                    // a1 = x0^2 + x1^2
    a1 = a1.fpSqrt ();                     // a1 = (x0^2 + x1^2)^(1/2)
    a1 = a1.fpSub (x0);                    // a1 = -x0 + (x0^2 + x1^2)^(1/2)

    if (a1.fpEquals (Felm.ZERO))
      return new F2elm (x0.fpSqrt (), Felm.ZERO);

    a1 = a1.fpDiv2 ();                     // a1 = [-x0 + (x0^2 + x1^2)^(1/2)] / 2
    a1 = a1.fpSqrt ();                     // a1 = ([-x0 + (x0^2 + x1^2)^(1/2)] / 2) ^ (1/2)
    
    a0 = a1.fpLeftShift (1);               // a0 = 2*a1
    a0 = a0.fpInverse ();                  // a0 = 1/(2*a1)
    a0 = x1.fpMult (a0);                   // a0 = x1 / (2*a1)
    
    return new F2elm (a0, a1);
  }
   

  private boolean isQuadraticResidue () {
    F2elm c;
    BigInteger d;

    d = Felm.getPrime ();
    d = d.multiply (d);
    d = d.subtract (BigInteger.ONE);
    d = d.shiftRight (1);                  // d = (p^2 - 1) / 2
    c = f2Pow (d);

    if (c.f2Equals (F2elm.ONE))
      return true;
    return false;
  }
    
    
  public F2elm f2Pow (BigInteger exp) {
    F2elm result = F2elm.ONE;
    F2elm base = new F2elm (this);
    BigInteger e = exp;

    while (e.compareTo (BigInteger.ZERO) == 1) {
      if (e.testBit (0) == true)
	result = result.f2Mult (base);
      e = e.shiftRight (1);
      base = base.f2Sqr ();
    }

    return result;
  }

    
  public F2elm f2Inverse () {
    // Compute inverse c = (x0 - x1*i) / (x0^2 + x1^2)
    Felm t0, t1, c0, c1;

    t0 = x0.fpSqr();                        // t0 = x0^2
    t1 = x1.fpSqr();                        // t1 = x1^2
    t0 = t0.fpAdd(t1);                      // t0 = x0^2 + x1^2
    t0 = t0.fpInverse();                    // t0 = 1 / (x0^2 + x1^2)
    t1 = x1.fpNegate();                     // t1 = -x1

    c0 = t0.fpMult(x0);                     // c0 = x0 / (x0^2 + x1^2)
    c1 = t1.fpMult(t0);                     // c1 = -x1 / (x0^2 + x1^2)

    return new F2elm (c0, c1);
  }


  public static F2elm[] inv4Way (F2elm z0, F2elm z1, F2elm z2, F2elm z3) {
    // Computes simultaneous inversion of 4 elements
    
    F2elm res[] = new F2elm[4];
    
    res[0] = z0.f2Mult(z1);                 // res0 = z0*z1
    res[1] = z2.f2Mult(z3);                 // res1 = z2*z3
    res[2] = res[0].f2Mult(res[1]);         // res2 = z0*z1*z2*z3
    res[3] = res[2].f2Inverse();            // res3 = 1/(z0*z1*z2*z3)
    res[2] = res[1].f2Mult(res[3]);         // res2 = 1/(z2*z3)
    res[3] = res[0].f2Mult(res[3]);         // res3 = 1/(z0*z1)
    res[0] = res[2].f2Mult(z1);             // res0 = 1/z0
    res[1] = res[2].f2Mult(z0);             // res1 = 1/z1
    res[2] = res[3].f2Mult(z3);             // res2 = 1/z2
    res[3] = res[3].f2Mult(z2);             // res3 = 1/z3

    return res;
  }


  public F2elm f2Swap (F2elm y, BigInteger option) {
    // Constant time swap regardless of whether option is 0 or 1
    Felm y0, y1;

    y0 = x0.fpSwap(y.x0, option);
    y1 = x1.fpSwap(y.x1, option);

    return new F2elm (y0, y1);
  }


  public F2elm f2Select (F2elm y, BigInteger option) {
    // Return this if option = 0 and y if option = 1
    BigInteger y0, y1, z0, z1, mask, bix0, bix1;

    mask = option.negate();                 // if option = 1 then mask = 1...1 because BigInteger
                                            // automatically sign extends as necessary

    // Get x0, x1, y0, y1 as their BigInteger representations
    y0 = (y.x0).fpGetValue();
    y1 = (y.x1).fpGetValue();
    bix0 = x0.fpGetValue();                 
    bix1 = x1.fpGetValue();                 

    z0 = bix0.xor(y0);                      // z0 = x0 xor y0
    z1 = bix1.xor(y1);                      // z1 = x1 xor y1

    z0 = z0.and(mask);                      // if mask = 0 then z0 = 0 else z0 = x0 xor y0 
    z1 = z1.and(mask);                      // if mask = 0 then z1 = 0 else z1 = x1 xor y1 

    z0 = bix0.xor(z0);                      // if mask = 0 then z0 = x0 
                                            // else z0 = x0 xor x0 xor y0 = y0

    z1 = bix1.xor(z1);                      // if mask = 0 then z1 = x1 
                                            // else z1 = x1 xor x1 xor y1 = y1

    return new F2elm (z0, z1);
  }


  public String toString() {
    return x1 + "*i + " + x0;
  }


  public byte[] toByteArray() {
    byte[] retval = new byte[2*Felm.primesize];
    System.arraycopy (x0.toByteArray(), 0, retval, 0, Felm.primesize);
    System.arraycopy (x1.toByteArray(), 0, retval, Felm.primesize, Felm.primesize);
    return retval;
  }
}


