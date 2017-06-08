
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

class InvalidFieldException extends Exception {
}


class Felm {
  /* Elements of GF(p) */

  // Set p to a placeholder value until the prime has been set
  public static BigInteger p = BigInteger.valueOf(2);
  private BigInteger value;
  public static int primesize = 1;    // number of bytes needed to represent the prime

  // Constants of the class
  public static final Felm ZERO = new Felm (BigInteger.ZERO);
  public static final Felm ONE = new Felm (BigInteger.ONE);


  public Felm (BigInteger v) {
    value = v.mod(p);
  }


  public Felm (Felm a) {
    value = a.fpGetValue();
  }


  public Felm (byte[] bytes) {
    value = new BigInteger (bytes);
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


  public Felm fpAdd (Felm y) {
    return new Felm (value.add (y.value));
  }


  public Felm fpSub (Felm y) {
    return fpAdd (y.fpNegate());             // this - y = this + (-y)
  }


  public Felm fpMult (Felm y) {
    return new Felm (value.multiply (y.value));
  }


  public Felm fpSqr() {
    return fpMult(this);
  }


  public boolean fpIsZero() {
    return value.equals(BigInteger.ZERO);
  }


  public boolean fpIsEven() {
    BigInteger c = value.and(BigInteger.ONE);
    return c.equals(BigInteger.ZERO);
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


  public Felm fpDiv2() {
    BigInteger v = value;

    if (fpIsOdd()) 
      v = v.add(p);                          // If this is odd, adding p makes it even
      v = v.shiftRight(1);                   // Divide by two = right shift by one

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
    x0 = a0;
    x1 = a1;
  }


  public F2elm (F2elm a) {
    x0 = a.x0;
    x1 = a.x1;
  }


  public F2elm (byte[] bytes) {
    int len = (bytes.length) / 2;
    x0 = new Felm (Arrays.copyOfRange (bytes, 0, len));
    x1 = new Felm (Arrays.copyOfRange (bytes, len, 2*len));
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


  public F2elm f2Sqr() {
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


  public boolean f2IsEven() {
    return x0.fpIsEven() && x1.fpIsEven();
  }


  public F2elm f2Div2() {
    Felm c0, c1;

    c0 = x0.fpDiv2();
    c1 = x1.fpDiv2();

    return new F2elm (c0, c1);
  }


  public F2elm f2Inverse() {
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


