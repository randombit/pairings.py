
This is a Python implementation of a pairing over a 256-bit BN curve.
It interops with the implementation described in
http://polycephaly.org/projects/dclxvi/ as well as the Golang x/crypto
bn256 module.

An elliptic curve pairing allows useful protocols such as identity-based
encryption, attribute-based encryption, group signatures, etc. For
details consult Eurocrypt proceedings and your pineal gland.

Warning: this curve was designed to provide 128-bit security, but due
to recent advances in the discrete logarithm in F_p^n it may provide
somewhat less than that. Estimates vary between 96 and 110 bits,

Extending this module to support a larger curve is an eventual goal; a
448-bit BN curve is currently estimated to be sufficient to provide
~128-bit security. Though alternative curve formulations such as BLS
may have superior performance at this security level.

Also, being Python this code is not particularly fast :) nor obviously
does it provide meaningful side channel protections. So you probably
do not want to use it for anything at all ever. It is primarily
intended for learning, and as a prototype for an implementation of the
same curve in C++ or Rust.

Enjoy!

Implemented:
 - [X] Optimal ate pairing over a 256-bit BN curve
 - [X] Point compression for G1
 - [X] Hashing to G1 (https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf)
 - [X] Example: BLS signature scheme
 - [X] Example: Tripartite key exchange

TODO:
 - [ ] Extend to support larger prime fields (eg BN-448)
 - [ ] Elligator hashing onto G1 (https://eprint.iacr.org/2014/043.pdf)
 - [ ] Point compression for G2
 - [ ] Example: Boneh-Franklin IBE

Some useful papers that I referenced in writing this code

"New software speed records for cryptographic pairings"
(https://cryptojedi.org/papers/dclxvi-20100714.pdf) describes the
curve. The line functions for the optimal ate pairing follow dclxvi.

"High-Speed Software Implementation of the Optimal Ate Pairing over
Barreto-Naehrig Curves" (https://eprint.iacr.org/2010/354) provided
most of the algorithms used in the field tower.

"Multiplication and Squaring on Pairing-Friendly Fields"
(https://eprint.iacr.org/2006/471)

"Pairing-Friendly Elliptic Curves of Prime Order"
(https://eprint.iacr.org/2005/133) is the paper introducing BN curves.

"Implementing Cryptographic Pairings over Barreto-Naehrig Curves"
(https://eprint.iacr.org/2007/390)

