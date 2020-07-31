import bn256
import random

# this is a literal implementation of the algorithms in Hierarchical Identity Based Encryption with Constant Size Ciphertext
# by Boneh, Boyen, and Goh
# https://crypto.stanford.edu/~dabo/papers/shibe.pdf

curve_order = bn256.order

# To generate system parameters for an HIBE of maximum depth l, select a random number generator g \belongs_to G, a random alpha \belongs_to Z_p,
# and set g_1 = g^alpha. Next, pick random elements g_2, g_3, h_1...h_l, \belongs_to G
# params = (g, g_1, g_2, g_3, h_1...h_l), where l = max_depth
# master_key = g_2^alpha
def Setup(l):
    k, g = bn256.g2_random()

    alpha = random.randrange(2, curve_order)

    # it turns out that gt.scalar_mul works for any kind of point... in other words g.scalar_mul(alpha) == bn256.curve_twist(g.x, g.y, g.z).scalar_mul(alpha)
    g_1 = g.scalar_mul(alpha)

    k, g_2 = bn256.g1_random()

    k, g_3 = bn256.g1_random()

    h_arr = list(map(lambda x: bn256.g1_random()[1], range(l)))

    master_key = g_2.scalar_mul(alpha)

    return (g, g_1, g_2, g_3, h_arr, master_key)

# To generate a private key d_{id} for identity ID = [I_1,...I_k] \belongs_to (Z_p^*)^k of depth k<=l, pick a random r \belongs_to Z_p and output
# d_{ID} =  () (master_key*(h_arr[0:k]))^r, g^r, h_arr[k+1]^r...h_arr[l]^r ) <--- 3-tuple
def KeyGen(params, id):
    (g, g_1, g_2, g_3, h_arr, master_key) = params

    depth = len(id)
    maxDepth = len(h_arr)

    r = random.randrange(2, curve_order)

    # compute first element
    running_product = bn256.curve_point(g_3.x, g_3.y, g_3.z)
    for x in range(0, depth):
        # h_x^{I_x}
        h_member_to_the_i = h_arr[x].scalar_mul(id[x])
        running_product = bn256.point_add(running_product, h_member_to_the_i)
    
    running_product = running_product.scalar_mul(r)

    result_1 = bn256.point_add(running_product, master_key)

    # compute second element
    result_2 = g.scalar_mul(r)

    # compute third element which is h_arr where every element is scalar_mul'd by r
    result_3 = list(map(lambda x: h_arr[x].scalar_mul(r), range(depth, maxDepth)))

    return (result_1, result_2, result_3)

# to encrypt a message M \belongs_to G_1 under the public key ID = (I_1, .. I_k) \belongs_to (Z_p^*)^k, pick a random s \belongs_to Z_p and output
# CT = (e(g_1, g_2)^s*m, g^s, (h_arr[0:k]*g_3)^s)
def Encrypt(params, id, m):
    (g, g_1, g_2, g_3, h_arr, master_key) = params

    depth = len(id)

    s = random.randrange(2, curve_order)

    # compute first element
    e = bn256.optimal_ate(g_1, g_2)
    e_to_s = e.exp(s)
    
    result_1 = e_to_s.mul(m)

    # compute second element
    result_2 = g.scalar_mul(s)

    # compute third element
    running_product = bn256.curve_point(g_3.x, g_3.y, g_3.z)
    for x in range(0, depth):
        # h_x^{I_x}
        h_member_to_the_i = h_arr[x].scalar_mul(id[x])
        running_product = bn256.point_add(running_product, h_member_to_the_i)

    running_product = running_product.scalar_mul(s)

    result_3 = running_product

    return (result_1, result_2, result_3)

# to decrypt a given ciphertext CT = (A, B, C) using the private key d_{ID} = (a_0, a_1, b_arr).
# output A*e(a_1,C) / e(B, a_0)
def Decrypt(private_key, ciphertext):
    (a_0, a_1, b_arr) = private_key
    (A, B, C) = ciphertext

    e_a1_C = bn256.optimal_ate(a_1, C)
    numerator = e_a1_C.mul(A)

    denominator = bn256.optimal_ate(B, a_0)

    return numerator.mul(denominator.inverse())

message = bn256.optimal_ate(bn256.g2_random()[1], bn256.g1_random()[1])

params = Setup(10)
path = [1,2,3]

# TODO: in order to vend keys for children these functions need to have params and master_key separated, so that children
# can generate further keys based on a_1
privkey = KeyGen(params, path)
ciphertext = Encrypt(params, path, message)
plaintext = Decrypt(privkey, ciphertext)

print("------------------------ MESSAGE:")
print(message)
print("------------------------ SETUP:")
print(params)
print("------------------------ PATH:")
print(path)
print("------------------------ PRIVATE KEY :")
print(privkey)
print("------------------------ CIPHERTEXT:")
print(ciphertext)
print("------------------------ PLAINTEXT:")
print(plaintext)

print("-------------------------------------------")
print("-------------------------------------------")
print("- Does MESSAGE equal PLAINTEXT ?? -")
print(message == plaintext)
