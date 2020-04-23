# Iteration and operation count for information-set decoding over F<sub>q</sub>

My paper [Information-set decoding for linear codes over
F<sub>q</sub>](http://eprint.iacr.org/2009/589) presents a new
algorithm for decoding linear codes over arbitrary finite fields
F<sub>q</sub>.


## Iteration count scripts

We use scripts to estimate the complexity of
information-set decoding attacks for *q*-ary codes and to
determine parameters for the McEliece Cryptosystem over
F<sub>q</sub>.

* [isdfq.c](isdfq.c) is written in C and uses the [MPFI library](https://gforge.inria.fr/projects/mpfi/).

* [isdfq.gp](isdfq.gp) is to be used with the
  [PARI/GP](https://pari.math.u-bordeaux.fr/) computer algebra system.
  The counts are a little less precise.

For details on how to work these scripts and a couple of examples,
see my website 
[https://cbcrypto.org/publications/scripts/](https://cbcrypto.org/publications/scripts/).

