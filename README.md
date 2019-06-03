This software is an implementation of the Supersingular Isogeny Diffie-Hellman (SIDH) key
exchange invented in 2011 by DeFeo, Jao, and Plut:

De Feo, Luca; Jao, David; Plut, Jerome. "Towards quantum-resistant cryptosystems from supersingular elliptic curve isogenies". PQCrypto 2011. Also available at https://eprint.iacr.org/2011/506.pdf.

This implementation is largely a Java translation of Microsoft's C implementation of SIDH except that the Java implementation is
object oriented. Their paper and code can be accessed using the following links:

1. https://eprint.iacr.org/2016/413
2. https://github.com/Microsoft/PQCrypto-SIDH

This implementation includes the parameter sets corresponding to the 434 bit prime (added in round 2 of the NIST competition), 
the 503 bit prime (currently the default option), and the 751 bit prime. Larger primes are more conservative from a security
perspective, but the smaller primes yield better performance.

In implementations of classical key exchanges, there is typically one function for key generation and one function for key 
agreement. This implementation requires separate versions for Alice and Bob because the two sides of the exchange are different.
If this code is used to replace a classical exchange in an existing application, be sure to call the correct functions for
key generation and key agreement (depending on whether they are being called by Alice or Bob).
