
package rlwe;

/**************************************************************************************************
 *
 * Implements elements over finite field GF(q) for use in RLWE implementation. The primality of q
 * is not checked because Fq need not be a valid field for the algorithm to work. In particular the 
 * BCNS implementation uses q = 2^32 - 1 which is not prime.
 *  
 * TODO: 
 *    - Write separate methods for faster computation of modular arithmetic for specific primes.
 *
 **************************************************************************************************/

import java.util.Arrays;
import java.lang.System;
import java.lang.Math;
import java.nio.ByteBuffer;

class Felm {

  // Set q to a placeholder value until the modulus has been set
  public static long q = 2;
  public static int numbits = 1;             // maximum number of bits in any Felm (depends on q)
  private long value;

  // Class constant
  public static final Felm ZERO = new Felm (0);
  public static final Felm ONE = new Felm (1);


  public Felm () {
    value = 0;
  }

  public Felm (long v) {
    // If v is negative v % q returns a negative number in Java.
    value = (v % q + q) % q;
  }


  public Felm (Felm a) {
    value = a.fqGetValue();
  }


  public static void setModulus (long qIn) {
    q = qIn;
    numbits = Long.SIZE - Long.numberOfLeadingZeros (q-1);
  }


  public static long getModulus () {
    return q;
  }


  public long fqGetValue () {
    return value;
  }


  public void fqSetValue (long v) {
    value = v;
  }


  public Felm fqAdd (Felm y) {
    return new Felm (value + y.value);
  }


  public Felm fqSub (Felm y) {
    return new Felm (value - y.value);
  }


  public Felm fqMult (Felm y) {
    return new Felm (value * y.value);
  }


  public Felm fqSqr() {
    return new Felm ((long) Math.pow(value, 2));
  }


  public boolean fqIsZero() {
    return value == 0;
  }


  public boolean fqIsEven() {
    return (value & 1) == 0;
  }


  public boolean fqIsOdd() {
    return (value & 1) == 1;
  }


  public boolean fqEqual (Felm y) {
    return value == y.value;
  }


  public boolean fqLessThan (Felm y) {
    return value < y.value;
  }


  public boolean fqGreaterThan (Felm y) {
    return value > y.value;
  }


  public Felm fqNegate() {
    return new Felm (q - value);
  }


  // Uses extended Euclidean algorithm
  public Felm fqInverse() {
    long a = value, m = q, x = 0, y = 1, quotient, tempy, tempa;

    while (a != 0) {
      quotient = m / a;
      tempy = y;
      y = x - quotient * y;
      x = tempy;
      tempa = a;
      a = m - quotient * a;
      m = tempa;
    }
    
    return new Felm (x);
  }


  public Felm fqDiv2() {
    long v = value;

    if (fqIsOdd()) 
      v = v + q;                             // If this is odd, adding q makes it even
    v = v >> 1;                              // Divide by two = right shift by one

    return new Felm (v);
  }


  public Felm fqDiv (Felm d) {
    Felm dinv = d.fqInverse();
    return fqMult (dinv);
  }


  public Felm rShift (int sb) {
    return new Felm (value >> sb);
  }


  public long getBit (int pos) {
    return value & (1 << pos);
  }


  public String toString() {
    return Long.toString(value);
  }


  public byte[] toByteArray() {
    return ByteBuffer.allocate(Long.SIZE).putLong(value).array(); 
  }
} 


