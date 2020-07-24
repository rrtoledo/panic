from bplib.bp import BpGroup
from gsproof import GSProof
import petlib.pack
from hashlib import sha256

from cca2_cs import CSEnc
from usps import USPS
from bsps import BSPS
from sok import SOK

class GSPO():
#Group Signature with Proof of Ownership
	def __init__(self):
		print("GSPO: Init")

class GSPO_auth():
	def __init__(self):
		print("GSPO_auth: Init")

	def GSetup(self):
		#print("GSPO_auth: GSetup")
		self.GS = GSProof()
		self.ck, self.xk, self.P = self.GS.ExtGen()
		self.order, self.G, self.u1, self.v1, self.w1, self.u2, self.v2, self.w2 = self.ck
		_, self.Xsi, self.Psi = self.xk
		self.order, self.G, self.g1, self.g2, self.e = self.P

		self.param_enc = (self.order, self.G, self.v1[0], self.v1[1])
		self.param_sig = (self.order, self.G, self.g1, self.g2, self.e)
		self.H = sha256

		self.Enc = CSEnc()
		self.ek, self.dk = self.Enc.key_gen(self.param_enc)

		self.Sig = BSPS()
		self.skI, self.vkI = self.Sig.keygen(self.param_sig)

		self.crs = (self.u1, self.v1, self.w1, self.u2, self.v2, self.w2)
		self.params = (self.P, self.crs, self.vkI, self.ek, self.H)
		return (self.params, self.skI, self.dk)


	def GRegister(self, vki):
		#print("GSPO_auth: GRegister")
		hi = self.g1 * self.order.random()
		sig0 = self.Sig.sign(self.param_sig, self.skI, [self.g1*0, hi], vki) 
		#print("GRegister, sig0")
		return hi, vki, sig0

	def GCertify(self, pki, t):
		#print("GSPO_auth: GCertify")
		hi, vki, sig0 = pki
		if self.Sig.verify(self.param_sig, self.vkI, [self.g1*0, hi], vki, sig0):
			sigt = self.Sig.sign(self.param_sig, self.skI, [self.g1*t, hi], vki)
			return t, sigt

class GSPO_prov():
	def __init__(self, params):
		self.P, self.crs, self.vkI, self.ek, self.H = params
		self.order, self.G, self.g1, self.g2, self.e = self.P
		self.param_sig = (self.order, self.G, self.g1, self.g2, self.e)
		self.param_sok = (self.order, self.G, self.g1)
		self.Sig = USPS()
		self.sok = SOK()

	def KeyGen(self):
		#print("GSPO_prov: KeyGen")
		self.sk, self.vk = self.Sig.keygen(self.param_sig)

	def GJoin(self, auth):
		#print("GSPO_prov: GJoin")
		self.hi, self.vk, self.sig0 = auth.GRegister(self.vk)
		self.pk = (self.hi, self.vk, self.sig0)

	def GJoint(self, auth, t):
		#print("GSPO_prov: GJoint")
		self.t, self.sigt = auth.GCertify(self.pk, t)

	def GIssue(self, m, pkui, pi_ui):
		#print("GSPO_prov: GIssue")
		if self.sok.verify(self.param_sok, pkui, self.hi, pi_ui):
			msg = [pkui]
			msg.extend(m)
			sig = self.Sig.sign(self.param_sig, self.sk, msg)
			return sig, (self.t, self.sigt)
		#print("--------------------------------------------- GSPO_prov: authentication error")
		return -1

	def GVerify(self, auth, pk_uv, pi_uv, pi_iuv):
		#print("GSPO_prov: GVerify")
		res = 0
		if self.sok.verify(self.param_sok,pk_uv, self.hi, pi_uv, pi_iuv):
			res = 1
			import time
			tverify = time.time()
			for i in range(len(pi_iuv)):
				#print("proof #"+str(i)+":", pi_iuv[i][1][0])
				proof = pi_iuv[i]
				pi, equation = proof
				pi2_v1, pi2_w1, pi1_v2, pi1_w2 = pi
				eq, X, Y, C, T = equation
				tt1 = time.time()
				boolean = auth.GS.Verify(eq, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
				#print("---", time.time()-tt1)
				if (boolean != 1):
					print("error proof", i, boolean)
				res *= boolean
			tverifyend = time.time()
			#print("proof verification time:", tverifyend - tverify)
		return res

class GSPO_user():
	def __init__(self, params):
		self.P, self.crs, self.vkI, self.ek, self.H = params
		self.order, self.G, self.g1, self.g2, self.e = self.P
		_, self.v1, _, _, _, _ = self.crs
		self.param_sig = (self.order, self.G, self.g1, self.g2, self.e)
		self.param_sok = (self.order, self.G, self.g1)
		self.param_enc = (self.order, self.G, self.v1[0], self.v1[1])
		self.Sig_auth = BSPS()
		self.Sig_prov = USPS()
		self.sok = SOK()
		self.CSEnc = CSEnc()

	def GKeyGen(self):
		#print("GSPO_user: GKeyGen")
		self.sku = self.order.random()

	def GDerive(self, pk):
		#print("GSPO_user: GDerive")
		h, vk, sig0 = pk
		if self.Sig_auth.verify(self.param_sig, self.vkI, [self.g1*0, h], vk, sig0):
			return h *self.sku

	def GAuthenticate(self, pk_x, pk_ux, m=""):
		#print("GSPO_user: Gauthenticate")
		h, vk, sig0 = pk_x

		if self.Sig_auth.verify(self.param_sig, self.vkI, [self.g1*0, h], vk, sig0) and pk_ux == h*self.sku:
			pi_ux = self.sok.prove(self.param_sok, pk_ux, h, self.sku, m)
			return pi_ux

	def GFinalize(self, auth, pki, m , sig_t, pkv, t):
		#print("GSPO_user: GFinalize")
		h_i, vk_i, sig0_i = pki
		h_v, vk_v, sig0_v = pkv
		sig_i, c_it = sig_t
		ti, sigt_i = c_it
		pk_ui = self.GDerive(pki)
		pk_uv = self.GDerive(pkv)

		bool_i0 = self.Sig_auth.verify(self.param_sig, self.vkI, [self.g1*0, h_i], vk_i, sig0_i)
		#print("bool_i0", bool_i0)

		bool_it = self.Sig_auth.verify(self.param_sig, self.vkI, [self.g1*t, h_i], vk_i, sigt_i)
		#print("bool_it", bool_it)

		bool_v = self.Sig_auth.verify(self.param_sig, self.vkI, [self.g1*0, h_v], vk_v, sig0_v)
		#print("bool_v", bool_v)

		msg_sig_i = [pk_ui]
		msg_sig_i.extend(m)
		bool_sig_i = self.Sig_prov.verify(self.param_sig, vk_i, msg_sig_i, sig_i)
		#print("bool_sig_i", bool_sig_i)

		if bool_i0 and bool_it and bool_v and bool_sig_i and ti==t:
			#print("--- all signatures verify")
			c_pk, r = self.CSEnc.enc(self.param_enc, auth.ek, h_i)

			from proofs_aggreg import prepare_proofs
			verify, pi_iuv = prepare_proofs(auth, self.vkI, pki, pkv, m, sig_t, t, pk_ui, pk_uv, self.sku, c_pk, auth.ek, r)

			#authentication proof on pk_uv
			pi_uv = self.GAuthenticate(pkv, pk_uv, pi_iuv)

			return pk_uv, pi_uv, pi_iuv


def main():
	print("--- init authority")
	auth = GSPO_auth()
	params, _, _ = auth.GSetup()

	print("\n--- init prover")
	prover = GSPO_prov(params)
	prover.KeyGen()
	prover.GJoin(auth)
	prover.GJoint(auth, 1)

	print("\n--- init verifier")
	verifier = GSPO_prov(params)
	verifier.KeyGen()
	verifier.GJoin(auth)

	print("\n--- init user")
	user = GSPO_user(params)
	user.GKeyGen()
	pk_ui = user.GDerive(prover.pk)
	pk_uv = user.GDerive(verifier.pk)

	print("\n--- user contacts prover")
	pi_ui = user.GAuthenticate(prover.pk, pk_ui)
	msg = [auth.G.hashG1(b"age=18")]
	sig_t = prover.GIssue(msg, pk_ui, pi_ui)

	print("\n--- user finalizes proof")
	pk_uv, pi_uv, pi_iuv = user.GFinalize(auth, prover.pk, msg, sig_t, verifier.pk, 1)

	print("\n--- user contacts verifier")
	boolean = verifier.GVerify(auth, pk_uv, pi_uv, pi_iuv)

	print("\n--- Access granted?", boolean)
	#TODO Trace

def time_benchmark(nb_test=100):
	print("--- Setup")
	auth = GSPO_auth()
	params, _, _ = auth.GSetup()

	prover = GSPO_prov(params)
	prover.KeyGen()
	prover.GJoin(auth)
	prover.GJoint(auth, 1)

	verifier = GSPO_prov(params)
	verifier.KeyGen()
	verifier.GJoin(auth)

	user = GSPO_user(params)
	user.GKeyGen()
	pk_ui = user.GDerive(prover.pk)
	pk_uv = user.GDerive(verifier.pk)

	pi_ui = user.GAuthenticate(prover.pk, pk_ui)
	msg = [auth.G.hashG1(b"age=18"), auth.G.hashG1(b"location=London")]
	sig_t = prover.GIssue(msg, pk_ui, pi_ui)

	time_final = []
	time_ver = []
	from time import time
	for i in range(nb_test):
		print(i, "/", nb_test)
		t_start = time()
		pk_uv, pi_uv, pi_iuv = user.GFinalize(auth, prover.pk, msg, sig_t, verifier.pk, 1)
		t_final = time()
		boolean = verifier.GVerify(auth, pk_uv, pi_uv, pi_iuv)
		t_verification = time()
		if boolean:
			time_final.append(t_final - t_start)
			time_ver.append(t_verification - t_final)
	print(len(time_final), "/" , n, "tests successful")
	import numpy as np
	final = np.array(time_final)
	ver = np.array(time_ver)
	print("Finalization:")
	print("- Average time:", np.mean(final))
	print("- Variance:", np.var(final))
	
	print("Verification:")
	print("- Average time:", np.mean(ver))
	print("- Variance:", np.var(ver))
