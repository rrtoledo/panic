from bplib.bp import BpGroup
from bplib.bp import G1Elem
from bplib.bp import G2Elem
from bplib.bp import GTElem

# Short Group Signatures via Structure-Preserving Signatures from Libert et al.
class USPS:

	def __init__(self, nb_msg = 6):
		print("SPS: Init")
		self.nb_msg = nb_msg

	def setup(self):
		"""Generates the Bilinear param"""
		#print("SPS: setup")
		G = BpGroup()
		g1, g2 = G.gen1(), G.gen2()
		e, o = G.pair, G.order()
		return (o, G, g1, g2, e)
	
	def keygen(self, params):
		"""
		Generates a fresh key pair
		private: {X_i}, {gamma_i}
		public: g_z, g_r, {g_i} = g_z^X_i g_r^gamma_i
		"""
		#print("SPS: KeyGen")
		(o, G, g1, g2, e) = params
		gz, gr = o.random()*g2, o.random()*g2
		sk = []
		pki=[]
		for i in range(self.nb_msg):
			xi, gammai=o.random(), o.random()
			sk.extend([[xi, gammai]])
			pki.extend([xi*gz+gammai*gr])
		pk = [gz, gr]
		pk.extend(pki)
	
		return (sk, pk)

	def sign(self, params, sk, m):
		"""Signs a list of group elements of G1"""
		#print("SPS: Sign")
		(o, G, g1, g2, e) = params
		if len(m)>self.nb_msg:
			print("Error: too many messages to sign")
			return (G1Elem.inf(G),G1Elem.inf(G))
		M = []
		for i in range(len(m)):
			if type(m[i]) == G1Elem:
				M.append(m[i])
			elif type(m[i]) == bytes:
				M.append(G.hashG1(m[i]))
			else:
				print("Error: type message")
				return -1
		sig = [sk[0][0].int_neg()*M[0], sk[0][1].int_neg()*M[0]]
		for i in range(1,len(M)):
			sig[0] = sig[0] + sk[i][0].int_neg()*M[i]
			sig[1] = sig[1] + sk[i][1].int_neg()*M[i]
	
		return sig

	def verify(self, params, pk, m, sig):
		"""
		Verifies a signature on messages m
		e(z, g_z) e(r, g_r) Pi e(M_i, g_i) == e(g1, g2)
		"""
		#print("SPS: Verify")
		(o, G, g1, g2, e) = params
		s0, s1 = sig
		if (s0 == G1Elem.inf(G) and s1== G1Elem.inf(G)):
			print("USPS: Verify --- Error: signature null")
			print("SPS: Verify", 0)
			return 0
		
		M = []
		for i in range(len(m)):
			if type(m[i]) == G1Elem:
				M.append(m[i])
			elif type(m[i]) == bytes:
				M.append(G.hashG1(m[i]))
			else:
				print("Error: type message")
				return -1
		
		ctr = 0
		for i in range(len(M)):
			if M[i] == G1Elem.inf(G):
				ctr+=1
		if ctr == len(M):
			print("USPS: Verify --- Error: message null")
			print("SPS: Verify", 0)
			return 0
		gz, gr = pk[0], pk[1]
		pki = pk[2:]
		res = e(s0,gz) * e(s1,gr)
		for i in range(len(M)):
			res = res * e(M[i],pki[i])
		return res.eq(GTElem.one(G))


	def proof_usps_hidesigandsigner(self, gsp, pk, M, sig):
		"""" creates GS proof that a USPS signature verifies
		with the verifying key, the signature and the first message secret"""
		#print("SPS: Prove")
		from bplib.bp import G1Elem
		from bplib.bp import G2Elem
		from bplib.bp import GTElem
		from petlib.bn import Bn
		params = gsp.P

		sps = USPS()
		gz, gr = pk[0], pk[1]
		pki = pk[2:]

		m = []
		for i in range(len(M)):
			if type(M[i]) == bytes:
				m.append(gsp.G.hashG1(M[i]))
			elif type(M[i]) == G1Elem:
				m.append(M[i])
			else:
				print("Error: wrong input, expected G1, got ", type(M[i]))
				return 0, []

		if len(m) < sps.nb_msg:
			for i in range(sps.nb_msg - len(m)):
				m.append(G1Elem.inf(gsp.G))

		if sps.verify(params, pk, m, sig) ==0:
			print("Signature does not verify")
			return 

		X = [{"type":"com", "value":sig[0]}, {"type":"com", "value":sig[1]}]
		B = [{"type":"pub", "value":G2Elem.inf(gsp.G)}, {"type":"pub", "value":G2Elem.inf(gsp.G)}]

		A = [{"type":"pub", "value":G1Elem.inf(gsp.G)}, {"type":"pub", "value":G1Elem.inf(gsp.G)}]
		Y = [{"type":"com", "value":gz}, {"type": "com", "value":gr}]

		for i in range(len(m)):
			if i == 0:
				X.append({"type":"pub", "value":m[i]})
			else:
				X.append({"type":"com", "value":m[i]})
			B.append( {"type":"pub", "value":G2Elem.inf(gsp.G)} )

		for j in range(len(pki)):
			Y.append({"type": "com", "value":pki[j]})
			A.append( {"type":"pub", "value":G1Elem.inf(gsp.G)} )

		C = []
		for i in range(len(X)):
			row = []
			for j in range(len(Y)):
				var = Bn(0)
				if i == j:
					var = Bn(1)
				row.append(var)
			C.append(row)

		success, res = gsp.CommitProof_eq("PPE", X, B, A, Y, C, GTElem.zero(gsp.G))
		verify = 0
		if success:
			eq_type, X1, Y1, C1, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res
			pi2_v1, pi2_w1, pi1_v2, pi1_w2 = gsp.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
			verify = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
			if verify:
				res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
		print("Do we successfully create a proof?", success)
		print("Does the proof successfully verify?", verify)
		print()

		return verify, res


# ---------- TESTS -------------

def test_setup():
	sps = USPS()
	sps.setup()

def test_keygen():
	sps = USPS()
	params = sps.setup()
	sk, pk = sps.keygen(params)

def test_sign():
	sps = USPS()
	params =  sps.setup()
	(o, G, g1, g2, e) = params
	sk, pk = sps.keygen(params)

	from petlib.bn import Bn
	from hashlib import sha256

	m = [G.hashG1(b"Hello World!"),G.hashG1(b"Hello me!"),G.hashG1(b"Hello us!")]
	sps.sign(params, sk,  m)

def test_verify():
	sps = USPS()
	params = sps.setup()
	(o, G, g1, g2, e) = params
	sk, pk = sps.keygen(params)

	from petlib.bn import Bn
	from hashlib import sha256

	m = [G.hashG1(b"Hello World!"),G.hashG1(b"Hello me!"),G.hashG1(b"Hello us!")]
	signature = sps.sign(params, sk,  m)

	assert sps.verify(params, pk, m, signature)

	m2 = [G.hashG1(b"Other Hello World!")]
	assert not sps.verify(params, pk, m2, signature)

def test_proof():
	from bplib.bp import G1Elem
	from bplib.bp import G2Elem
	from bplib.bp import GTElem
	from petlib.bn import Bn
	from gsproof import GSProof
	gsp = GSProof()
	gsp.ExtGen()
	params = gsp.P

	sps = USPS()
	sk, pk = sps.keygen(params)
	gz, gr = pk[0], pk[1]
	pki = pk[2:]

	M = [b"Hello World!", b"Hello me!", b"Hello us!"]
	m = []
	for i in range(len(M)):
		if type(M[i]) == bytes:
			m.append(gsp.G.hashG1(M[i]))
		elif type(M[i]) == G1Elem:
			m.append(M[i])
		else:
			print("Error: wrong input, expected G1, got ", type(M[i]))
			return 0, []

	if len(m) < sps.nb_msg:
		for i in range(sps.nb_msg - len(m)):
			m.append(G1Elem.inf(gsp.G))
	elif sps.nb_msg < len(m):
		return
	sig = sps.sign(params, sk,  m)

	verify, res = sps.proof_usps_hidesigandsigner(gsp, pk, M, sig)
	assert verify

