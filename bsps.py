from bplib.bp import BpGroup
from bplib.bp import G1Elem
from bplib.bp import G2Elem
from bplib.bp import GTElem
from petlib.bn import Bn


# Optimal Structure-Preserving Signature in Assymetric Bilinear Group
# Abe Groth et al.

class BSPS:
	msg_g1 = 2 #nb of messages in G1
	msg_g2 = 8 #nb of messages in G2

	def __init__(self):
		print("BSPS: Init")

	def setup(self):
		"""Generates the Bilinear param"""
		#print("BSPS: setup")
		G = BpGroup()
		g1, g2 = G.gen1(), G.gen2()
		e, o = G.pair, G.order()
		return (o, G, g1, g2, e)
	
	def keygen(self, params):
		"""Generates a fresh key pair"""
		#print("BSPS: KeyGen")
		o, _, g1, g2, _ = params
		v, z = o.random(), o.random()
		V = v*g2
		Z = z*g2
		u=[]
		U=[]
		w=[]
		W=[]

		for i in range(self.msg_g2):
			r=o.random()
			u.extend([r])
			U.extend([r*g1])
		for j in range(self.msg_g1):
			s=o.random()
			w.extend([s])
			W.extend([s*g2])
		sk = [u,w,v,z]
		pk = [U,W,V,Z]
	
		return (sk, pk)

	def sign(self, params, sk, m, n):
		"""Signs a list of group elements of G1 and G2"""
		#print("BSPS: Sign")
		o, G, g1, g2, _ = params
		if len(m)>self.msg_g1 or len(n)>self.msg_g2:
			print("BSPS: Sign --- Error: message(s) too long", len(m),
				">?", self.msg_g1, len(n), ">?", self.msg_g2)
			return (G1Elem.inf(G),G1Elem.inf(G),G2Elem.inf(G))
		u,w,v,z = sk

		r = o.random()
		R = r * g1
		temp = z.mod_sub(r.mod_mul(v,o),o)
		S = temp*g1
		T = g2
		
		for i in range(self.msg_g1):
			if i< len(m):
				S=S+w[i].int_neg()*m[i]
			else:
				S=S+w[i].int_neg()*g1
		for j in range(self.msg_g2):
			if j< len(n):
				T=T+u[j].int_neg()*n[j]
			else:
				T=T+u[j].int_neg()*g2

		return [R, S, r.mod_inverse(o)*T]

	def verify(self, params, pk, m, n, sig):
		"""
		Verifies a signature on messages m
		e(R,V) e(S,H) Pi e(M_i,W_i) == e(G,Z)
		e(R,T) Pi e(U_i, N_i) == e(G,H)
		"""
		#print("BSPS: Verify")
		R, S, T = sig
		o, G, g1, g2, e = params
		if (R.eq(G1Elem.inf(G)) and S.eq(G1Elem.inf(G)) and T.eq(G2Elem.inf(G))):
			print("BSPS: Verify --- Error: signature null")
			print("BSPS: Verify", 0)
			return 0

		U, W, V, Z = pk

		res1 = e(R,V) * e(S,g2)
		for i in range(self.msg_g1):
			if i< len(m):
				res1 = res1 * e(m[i],W[i])
			else:
				res1 = res1 * e(g1,W[i])

		res2 = e(R,T)
		for j in range(self.msg_g2):
			if j< len(n):
				res2 = res2 * e(U[j],n[j])
			else:
				res2 = res2 * e(U[j],g2)

		return res1.eq(e(g1,Z)) and res2.eq(e(g1,g2))

	
# ---------- TESTS -------------

def test_setup():
	sps = BSPS()
	sps.setup()

def test_keygen():
	sps = BSPS()
	params = sps.setup()
	sk, pk = sps.keygen(params)

def test_sign():
	sps = BSPS()
	params = sps.setup()
	o, G, g1, g2, e = params
	sk, pk = sps.keygen(params)

	from petlib.bn import Bn

	m = [G.hashG1(b"Hello World!"), G.hashG1(b"Hello me!")]
	n = [Bn.from_binary(b"Hello you!").mod_add(Bn(0),o) * g2]
	sig = sps.sign(params, sk,  m, n)
	print(sig)

def test_verify():
	sps = BSPS()
	params = sps.setup()
	o, G, g1, g2, e = params
	sk, pk = sps.keygen(params)

	from petlib.bn import Bn

	m = [G.hashG1(b"Hello World!"),G.hashG1(b"Hello me!")]
	n = [g2]
	signature = sps.sign(params, sk, m, n)
	print("verify true")
	assert sps.verify(params, pk, m, n, signature)

	m2 = [G.hashG1(b"Other Hello World!")]
	n2 = [Bn.from_binary(b"Other Hello you!").mod_add(Bn(0),o) * g2]
	print("verify false")
	assert not sps.verify(params, pk, m2, n2, signature)

