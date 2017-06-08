package rlwe;

/**************************************************************************************************
 *
 * Implements a class for polynomials with Felm coefficients and a class for elements of a 
 * polynomial quotient ring.
 *
 * TODO:
 *    - Implement faster polynomial multiplication based on FFT.
 * 
 **************************************************************************************************/

import java.math.BigInteger;
import java.util.Arrays;
import java.lang.System;
import java.lang.Math;
import java.nio.ByteBuffer;
import java.nio.LongBuffer;


class Polynomial {

  public int degree;
  public Felm[] coeff;


  public Polynomial () {
    degree = 0;
    coeff = null;
  }

  public Polynomial (int d) {
    degree = d;
    coeff = new Felm[degree];

    for (int i = 0; i < degree; i++)
      coeff[i] = new Felm();
  }


  public Polynomial (long[] coefficient) {
    degree = coefficient.length;
    coeff = new Felm[degree];
    for (int i = 0; i < degree; i++)
      coeff[i] = new Felm (coefficient[i]);
  } 


  public Polynomial (Felm[] coefficient) {
    degree = coefficient.length;
    coeff = new Felm[degree];
    for (int i = 0; i < degree; i++) 
      coeff[i] = new Felm (coefficient[i]);
  }


  public Polynomial (Polynomial b) {
    degree = b.degree;
    coeff = new Felm[degree];

    for (int i = 0; i < degree; i++)
      coeff[i] = new Felm (b.coeff[i]);    
  }


  public Polynomial (byte[] inBytes) {
    LongBuffer lb = ByteBuffer.allocate (inBytes.length).put(inBytes).asLongBuffer();
    degree = lb.limit();
    for (int i = 0; i < degree; i++)
      coeff[i] = new Felm (lb.get (i));
  }


  public int getDegree() {
    return degree;
  }
  
   
  public Felm getCoeff (int index) {
    if (index >= degree || index < 0)
      return Felm.ZERO;
    return coeff[index];
  }


  public void setCoeff (int index, Felm value) {
    if (degree > index && index >= 0)
      coeff[index] = value;
  }


  public void setCoeff (int index, long value) {
    setCoeff (index, new Felm(value));
  }


  public boolean equal (Polynomial p) {
    for (int i = 0; i < Math.max (degree, p.degree); i++) 
      if (getCoeff (i) != p.getCoeff (i))
        return false;
    return true;
  }


  public Polynomial polyAdd (Polynomial a) {
    Polynomial c = new Polynomial (Math.max (a.degree, degree));
    int adiff = c.degree - a.degree;
    int tdiff = c.degree - degree;

    for (int i = 0; i < c.degree; i++) {
      Felm ci = getCoeff(i-tdiff).fqAdd (a.getCoeff(i-adiff));
      c.setCoeff (i, ci);
    }

    return c;
  }


  public Polynomial polySub (Polynomial a) {
    Polynomial c = new Polynomial (Math.max (a.degree, degree));

    for (int i = 0; i < c.degree; i++) {
      Felm ci = getCoeff(i).fqSub (a.getCoeff(i));
      c.setCoeff (i, ci);
    }

    return c;
  }


  public Polynomial polyMult (Polynomial a) {
    // polyMult calls either polySlowMult or polyFFTMult. 
    return polySlowMult (a);
  }


  public Polynomial polySlowMult (Polynomial a) {
    // Slow polynomial multiplication. 

    int i, j, cd = a.degree + degree;
    Felm singlepair;

    Felm[] c = new Felm[a.degree + degree - 1];
    for (i = 0; i < c.length; i++)
      c[i] = Felm.ZERO;    

    for (i = 0; i < a.degree; i++)
      for (j = 0; j < degree; j++) {
        singlepair = getCoeff (i).fqMult (a.getCoeff(j));
        c[i+j] = c[i+j].fqAdd (singlepair);
      }

    return new Polynomial (c);
  }


  public Polynomial polyFFTMult (Polynomial a) {
    // Fast polynomial multiplication based on FFT ... currently unimplemented
    return new Polynomial ();
  }


  public Polynomial polyMultAdd (Polynomial a, Polynomial b) {
    // Multiply by a and add b
    return b.polyAdd (this.polyMult (a));
  }


  public Polynomial polyMod (Polynomial m) {
    Polynomial c = new Polynomial (coeff);
    Felm multiplier, lead = m.getCoeff(0), newval;
    int sizediff = degree - m.degree;        

    for (int i = 0; i <= sizediff; i++) {
      multiplier =  c.getCoeff(i).fqDiv (lead);
      for (int j = i; j < m.degree+i; j++) {
        newval = multiplier.fqMult (m.getCoeff(j-i));
        newval = c.getCoeff (j).fqSub (newval);
        c.setCoeff (j, newval);
      }
    }

    return new Polynomial (Arrays.copyOfRange (c.coeff, sizediff+1, degree));
  }


  public Polynomial pointwiseMult (Polynomial a) {
    Polynomial c = new Polynomial (degree);
    for (int i = 0; i < degree; i++)
      c.setCoeff (i, coeff[i].fqMult (a.getCoeff(i)));
    return c;
  }


  public Polynomial pointwiseMultAdd (Polynomial a, Polynomial b) {
    Polynomial c = new Polynomial (degree);
    Felm t;

    for (int i = 0; i < degree; i++) {
      t = coeff[i].fqMult (a.getCoeff(i));
      t = t.fqAdd (b.getCoeff(i));
      c.setCoeff (i, t);
    }

    return c;
  }


  public Polynomial fftForward () {
    // fill this in
    return new Polynomial();
  }


  public String toString () {
    String s = "(" + getCoeff(0);

    for (int i = 1; i < degree; i++)
      s = s + ", " + getCoeff(i);

    return s + ")";
  }


  public byte[] toByteArray () {
    ByteBuffer bb = ByteBuffer.allocate (degree * Long.SIZE);
    for (int i = 0; i < degree; i++) 
      bb.putLong (i, coeff[i].fqGetValue());
    return bb.array();
  }
}



class RingElt extends Polynomial {
  // Implement elements of Fq[x] / m(x), elements of a polynomial ring mod a polynomial

  public static Polynomial modulus;
  public static boolean modSet = false;   // Indicates whether the modulus has been set


  public RingElt () {
    super (modSet? modulus.getDegree() - 1: 0);
  }

  
  public RingElt (int d) {
    super (modSet? d: 0);
  }


  public RingElt (long[] coefficient) {
    super (modSet? coefficient: null);
    reduce();
  } 


  public RingElt (Felm[] coefficient) {
    super (modSet? coefficient: null);
    reduce();
  }


  public RingElt (Polynomial b) {
    super (modSet? b: null);
    reduce();
  }


  public RingElt (RingElt r) {
    super (modSet? r.coeff : null);
  }


  public RingElt (byte[] inBytes) {
    super (modSet? inBytes: null);
    reduce();
  }


  public static Polynomial getModulus () {
    return modulus;
  }


  public static void setModulus (Polynomial q) {
    modulus = new Polynomial (q);
    modSet = true;
  }


  public RingElt ringAdd (RingElt r) {
    return new RingElt (polyAdd (r).polyMod (modulus));
  }


  public RingElt ringSub (RingElt r) {
    return new RingElt (polySub (r).polyMod (modulus));
  }


  public RingElt ringMult (RingElt r) {
    return new RingElt (polyMult (r).polyMod (modulus));
  }


  public RingElt ringMultAdd (RingElt r1, RingElt r2) {
    return r2.ringAdd (r1.ringMult(this));
  }


  private void reduce() {
    Polynomial p = new Polynomial (coeff);
    coeff = p.polyMod (modulus).coeff;
    degree = coeff.length;
  }
}
