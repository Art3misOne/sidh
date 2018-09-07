This software is an implementation of the Supersingular Isogeny Diffie-Hellman (SIDH) key
exchange invented in 2011 by DeFeo, Jao, and Plut:

De Feo, Luca; Jao, David; Plut, Jerome. "Towards quantum-resistant cryptosystems from supersingular elliptic curve isogenies". PQCrypto 2011. Also available at https://eprint.iacr.org/2011/506.pdf.

This implementation is largely a Java translation of Microsoft's C implementation of SIDH except that the Java implementation is
object oriented. Their paper and code can be accessed using the following links:

https://eprint.iacr.org/2016/413

https://github.com/Microsoft/PQCrypto-SIDH

The latest update to this code incorporates the optimizations in version 3 of Microsoft's code and includes the additional 
parameter set using a 503 bit prime.

In implementations of classical key exchanges, there is typically one function for key generation and one function for key 
agreement. This implementation requires separate versions for Alice and Bob because the two sides of the exchange are different.
If this code is used to replace a classical exchange in an existing application, be sure to call the correct functions for
key generation and key agreement (depending on whether they are being called by Alice or Bob).
