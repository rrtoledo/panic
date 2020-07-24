from bplib.bp import BpGroup
from hashlib import sha256
from petlib.bn import Bn
import pytest

class CCA2EG:

	def __init__(self):
		print("CCA2EG: Init")

	def params_gen(self):
		"""Generates the AHEG for an EC group nid"""
		print("CCA2EG: Setup")
		G = BpGroup()
		g1, g2 = G.gen1(), G.gen2()
		e, o = G.pair, G.order()
		return (G, g1, o)

	def key_gen(self, params):
		"""Generates a fresh key pair"""
		print("CCA2EG: Key Gen")
		_, g1, o = params
		priv1 = o.random()
		priv2 = o.random()
		priv = (priv1, priv2)
		pub1 = priv1 * g1
		pub2 = priv2 * g1
		pub = (pub1, pub2)
		return (pub, priv)

	def challenge(self, els):
		"""Packages a challenge in a bijective way"""
		elem = [len(els)] + els
		elem_str = map(str, els)
		elem_len = map(lambda x: "%s||%s" % (len(x) , x), elem_str)
		state = "|".join(elem_len)
		H = sha256()
		H.update(state.encode("utf8"))
		return H.digest()
    

	def enc(self, params, pub, m):
		"""Encrypts the values of a group element"""
		print("CCA2EG: Enc")
		G, g1, o = params

		r1 = o.random()
		r2 = o.random()
		c11 = r1 * pub[0] + m
		c21 = r2 * pub[1] + m
		c12 = r1 * g1
		c22 = r2 * g1

		d = o.random() * g1
		s1 = o.random()
		s2 = o.random()
		e11 = s1 * pub[0] + d
		e21 = s2 * pub[1] + d
		e12 = s1 * g1
		e22 = s2 * g1

		state = [g1, pub[0], pub[1], c11, c21, c12, c22, e11, e21, e12, e22]
		hash_c = self.challenge(state)
		c = Bn.from_binary(hash_c) % o
		z = c*m +d
		z1 = (r1 * c + s1) % o
		z2 = (r2 * c + s2) % o

		return ((c11, c12), (c21, c22), (e11, e12), (e21, e22), c, z, z1, z2)


	def dec(self, params, priv, pub, cipher):
		"""Decrypts an encrypted group element"""
		print("CCA2EG: Make")
		_, g1, o = params
		(c11, c12), (c21, c22), (e11, e12), (e21, e22), c, z, z1, z2 = cipher

		b = 1
		if c * c11 + e11 != z1 * pub[0] + z or  c * c12 + e12 != z1 * g1:
			b = 0
		if c * c21 + e21 != z2 * pub[1] + z or  c * c22 + e22 != z2 * g1:
			b = 0

		plain1 = c11 + (-priv[0] * c12)
		plain2 = c21 + (-priv[1] * c22)

		return plain1, plain2, b


def test_CCA2EG():
	eg = CCA2EG()
	params = eg.params_gen()
	G, g1, o = params
	(pub, priv) = eg.key_gen(params)

	# Check encryption and decryption
	m = o.random() * g1
	c1 = eg.enc(params, pub, m)
	m1, m2, b = eg.dec(params, priv, pub, c1) 
	print (m1, m2, m1==m2, b)
	assert b == 1 and m1 == m

