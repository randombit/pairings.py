#!/usr/bin/python

import random
import bn256

a_k = random.randrange(2, bn256.order)
b_k = random.randrange(2, bn256.order)
c_k = random.randrange(2, bn256.order)

a_g1 = bn256.g1_scalar_base_mult(a_k)
b_g1 = bn256.g1_scalar_base_mult(b_k)
c_g1 = bn256.g1_scalar_base_mult(c_k)

a_g2 = bn256.g2_scalar_base_mult(a_k)
b_g2 = bn256.g2_scalar_base_mult(b_k)
c_g2 = bn256.g2_scalar_base_mult(c_k)

a_key = bn256.gt_scalar_mult(bn256.optimal_ate(b_g2, c_g1), a_k)
b_key = bn256.gt_scalar_mult(bn256.optimal_ate(c_g2, a_g1), b_k)
c_key = bn256.gt_scalar_mult(bn256.optimal_ate(a_g2, b_g1), c_k)

print(bn256.gt_hash(a_key))
print(bn256.gt_hash(b_key))
print(bn256.gt_hash(c_key))
