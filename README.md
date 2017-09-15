This software is an implementation of the Supersingular Isogeny Diffie-Hellman (SIDH) key
exchange invented in 2011 by DeFeo, Jao, and Plut. It is largely a Java translation of
Microsoft's C language implementation of SIDH except that the Java implementation is
object oriented. This code is designed to be compatible for use in the Signal Messaging
app.

Note that in the Signal API there is one function to calculate the shared agreement and
one function to generate a key pair regardless of whether these functions are invoked by
the responder or initiator. These are different on each side in SIDH because the degrees of
the isogenies are different. Thus, this implementation has separate functions for computing
the shared agreement for the initiator and responder, and separate functions for generating
key pairs. These differences, though minor, need to be accounted for when integrating this
implementation into any existing application.
