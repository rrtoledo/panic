from __future__ import print_function
import time

from bplib.bp import BpGroup
from bplib.bp import G1Elem
from bplib.bp import G2Elem
from bplib.bp import GTElem
from petlib.bn import Bn

from usps import USPS
from bsps import BSPS


class GSProof():
	cms = ["bas","pub","enc","com","unt","sca"]

	eqs = {
		"PPE":{1:["bas","pub","enc","com"],2:["bas","pub","enc","com"]},
		"PN1":{1:["bas","pub","enc"],2:["bas","com"]},
		"PC1":{1:["bas","pub"],2:["bas","com"]},
		"PN2":{1:["bas","com"],2:["bas","pub","enc"]},
		"PC2":{1:["bas","com"],2:["bas","pub"]},
		"ME1":{1:["bas","pub","enc","com"],2:["unt", "sca"]},
		"MN1":{1:["bas","pub","enc"],2:["unt", "sca"]},
		"MC1":{1:["bas","pub"],2:["unt", "sca"]},
		"ML1":{1:["bas","com"],2:["unt"]},
		"ME2":{1:["unt", "sca"],2:["bas","pub","enc","com"]},
		"MN2":{1:["unt", "sca"],2:["bas","pub","enc"]},
		"MC2":{1:["unt", "sca"],2:["bas","pub"]},
		"ML2":{1:["unt"],2:["bas","com"]},
		"QE":{1:["unt", "sca"],2:["unt", "sca"]},
		"QC1":{1:["unt"],2:["unt", "sca"]},
		"QC2":{1:["unt", "sca"],2:["unt"]}
	}

	def __init__(self):
		#print("GS: Init")
		return
		

	def ExtGen(self):
		#print("GS: Setup Binding")
		self.G = BpGroup()
		self.e = self.G.pair
		self.order = self.G.order()
		self.g1, self.g2 = self.G.gen1(), self.G.gen2()
		self.gt = self.e(self.g1, self.g2)
		self.P = (self.order, self.G, self.g1, self.g2, self.e)

		rho = self.order.random()
		xsi = self.order.random()
		self.v1 =  [xsi * self.g1, self.g1]
		self.w1 =  [rho * xsi * self.g1, rho * self.g1]
		self.u1 =  [rho * xsi * self.g1, rho * self.g1 + self.g1]
		self.Xsi = [xsi.mod_inverse(self.order).int_neg(),1]

		sigma = self.order.random()
		psi = self.order.random()
		self.v2 =  [psi * self.g2, self.g2]
		self.w2 =  [sigma * xsi * self.g2, sigma * self.g2]
		self.u2 =  [sigma * xsi * self.g2, sigma * self.g2 + self.g2]
		self.Psi = [psi.mod_inverse(self.order).int_neg(),1]

		self.ck = (self.order, self.G, self.u1, self.v1, self.w1, self.u2, self.v2, self.w2)
		self.xk = (self.ck, self.Xsi, self.Psi)
		
		return (self.ck, self.xk, self.P)

	def G1pub(self, x1):
		val = {}
		val["type"] = "pub"
		val["group"] = 1
		val = self.CommitG1Group(val, x1, 0, 0)
		return val

	def G1enc(self, x1):
		r = self.order.random()
		val = {}
		val["type"] = "enc"
		val["group"] = 1
		val = self.CommitG1Group(val, x1, r, 0)
		return val

	def G1com(self, x1):
		r = self.order.random()
		s = self.order.random()
		val = {}
		val["type"] = "com"
		val["group"] = 1
		val =  self.CommitG1Group(val, x1, r, s)
		return val

	def G1bas(self):
		val = {}
		val["type"] = "bas"
		val["group"] = 1
		val = self.CommitG1Group(val, self.g1, 0, 0)
		return val

	def G1sca(self, x):
		r = self.order.random()
		val = {}
		val["type"] = "sca"
		val["group"] = 1
		val = self.CommitG1Scalar(val, x, r)
		return val

	def G1unt(self):
		val = {}
		val["type"] = "unt"
		val["group"] = 1
		val = self.CommitG1Scalar(val, 1, 0)
		return val

	def CommitG1Group(self, val, x1, r, s):
		c0 =      r*self.v1[0] + s*self.w1[0]
		c1 = x1 + r*self.v1[1] + s*self.w1[1]
		val["value"] = [c0, c1]
		val["committed"] = {1: x1, 2: G2Elem.inf(self.G), "Zp": 0, "r": r, "s": s}
		return val

	def CommitG1Scalar(self, val, x, r):
		c0 = x*self.u1[0] + r*self.v1[0]
		c1 = x*self.u1[1] + r*self.v1[1]
		val["value"] = [c0, c1]
		val["committed"] = {1: G1Elem.inf(self.G), 2: G2Elem.inf(self.G), "Zp": x, "r": r, "s":0}
		return val

	def CommitG1(self, x):
		try:
			ttype = x["type"]
			if ttype == "unt":
				return self.G1unt()
			elif ttype == "bas":
				return self.G1bas()

			value = x["value"]
			
			if type(value) == G1Elem:
				if ttype == "pub":
					return self.G1pub(value)
				elif ttype == "enc":
					return self.G1enc(value)
				elif ttype == "com":
					return self.G1com(value)
			
			if type(value) == Bn or type(value) == int:
				if type(value) == int:
					value = Bn(value)
				if ttype == "sca":
					return self.G1sca(value)

		except Exception as e:
			print("Error G1 commit", e)

	def G2pub(self, x2):
		val = {}
		val["type"] = "pub"
		val["group"] = 2
		val = self.CommitG2Group(val, x2, 0, 0)
		return val

	def G2enc(self, x2):
		r = self.order.random()
		val = {}
		val["type"] = "enc"
		val["group"] = 2
		val = self.CommitG2Group(val, x2, r, 0)
		return val

	def G2com(self, x2):
		r = self.order.random()
		s = self.order.random()
		val = {}
		val["type"] = "com"
		val["group"] = 2
		val = self.CommitG2Group(val, x2, r, s)
		return val

	def G2bas(self):
		val = {}
		val["type"] = "bas"
		val["group"] = 2
		val = self.CommitG2Group(val, self.g2, 0, 0)
		return val

	def G2sca(self, x):
		r = self.order.random()
		val = {}
		val["type"] = "sca"
		val["group"] = 2
		val = self.CommitG2Scalar(val, x, r)
		return val

	def G2unt(self):
		val = {}
		val["type"] = "unt"
		val["group"] = 2
		val = self.CommitG2Scalar(val, 1, 0)
		return val

	def CommitG2Group(self, val, x2, r, s):
		c0 =      r*self.v2[0] + s*self.w2[0]
		c1 = x2 + r*self.v2[1] + s*self.w2[1]
		val["value"] = [c0, c1]
		val["committed"] = { 1: G1Elem.inf(self.G), 2: x2, "Zp": 0, "r": r, "s": s }
		return val

	def CommitG2Scalar(self, val, x, r):
		c0 = x*self.u2[0] + r*self.v2[0]
		c1 = x*self.u2[1] + r*self.v2[1]
		val["value"] = [c0, c1]
		val["committed"] = { 1: G1Elem.inf(self.G), 2: G2Elem.inf(self.G), "Zp" : x, "r" : r, "s": 0}
		return val

	def CommitG2(self, y):
		try:
			ttype = y["type"]

			if ttype == "bas":
				return self.G2bas()
			elif ttype == "unt":
				return self.G2unt()
			
			value = y["value"]

			if type(value) == G2Elem:
				if ttype == "pub":
					return self.G2pub(value)
				elif ttype == "enc":
					return self.G2enc(value)
				elif ttype == "com":
					return self.G2com(value)
			
			if type(value) == Bn or type(value) == int:
				if type(value) == int:
					value = Bn(value)
				if ttype == "sca":
					return self.G2sca(value)
					
		except Exception as e:
			print("Error G2 commit", e)

	def Commit(self, to_commit):
		try:
			group = to_commit["group"]
			if group == 1:
				return self.CommitG1(to_commit)
			elif group == 2:
				return self.CommitG2(to_commit)
			else:
				print("Commit error: group to commit with unknown")
		except Exception as e:
			print("Error commit:", e)

	def CheckEqFormat(self, eq, X, Y, C):
		#print("CheckEqFormat")
		m = len(X)
		n = len(Y)

		if len(C) != m:
			print("Wrong size (#rows)")
			return 0
		for i in range(m):
			if len(C[i]) != n:
				print("Wrong size (#cols)")
				return 0

		potential_x_format_error = 0
		if m != 0 :
			one = list(set(x["type"] for x in X))
			for tone in one:
				if tone not in self.eqs[eq][1]:
					potential_x_format_error = 1
		if potential_x_format_error != 0:
			for i in range(m):
				if X[i]["type"] not in self.eqs[eq][1]:
					for j in range(n):
						if C[i][j] != 0:
							print("X wrong format:", tone, "not in", self.eqs[eq][1])
							return 0

		potential_y_format_error = 0
		if n != 0:
			two = list(set(y["type"] for y in Y))
			for ttwo in two:
				if ttwo not in self.eqs[eq][2]:
					potential_y_format_error = 1
		if potential_y_format_error != 0:
			for j in range(n):
				if Y[j]["type"] not in self.eqs[eq][2]:
					for i in range(m):
						if C[i][j] != 0:
							print("Y wrong format:", ttwo, "not in", self.eqs[eq][2])
							return 0

		if eq == "PPE":
			for i in range(m):
				x=X[i]
				if x["type"] == "pub" or x["type"] == "enc":
					for j in range(n):
						y = Y[j]
						if y["type"] == "pub" or y["type"] == "enc":
							if C[i][j]!= 0:
								print("cij = 0 for pub/enc vars", i,j)
								return 0
					
		return 1


	def CheckEqType(self, eq, X, Y, C):
		#print("CheckEqType")
		if eq in ["PPE", "PN1", "PC1", "PN2", "PC2"]:
			for x in X:
				if type(x["value"][0]) != G1Elem:
					print("element in C not from G1", x)
					return 0
				if "committed" in x:
					if x["committed"][2] != G2Elem.inf(self.G) or x["committed"]["Zp"] != 0:
						print("not a commitment of G1")
						return 0
			for y in Y:
				if type(y["value"][0]) != G2Elem:
					print("element in D not from G2")
					return 0
				if "committed" in y:
					if y["committed"][1] != G1Elem.inf(self.G) or y["committed"]["Zp"] != 0:
						print("not a commitment of G2")
						return 0

		elif eq in ["ME1", "MN1", "MC1", "ML1"]:
			for x in X:
				if type(x["value"][0]) != G1Elem:
					print("element in C not from G1", type(x["value"][0]))
					return 0
				if "committed" in x:
					if x["committed"][2] != G2Elem.inf(self.G) or x["committed"]["Zp"] != 0:
						print("not a commitment of G1")
						return 0
			for y in Y:
				if type(y["value"][0]) != G2Elem:
					print("element in D not from Zp", type(y["value"][0]))
					return 0
				if "committed" in y:
					if y["committed"][1] != G1Elem.inf(self.G) or y["committed"][2] != G2Elem.inf(self.G):
						print("not a commitment of Zp")
						return 0
		elif eq in ["ME2", "MN2", "MC2", "ML2"]:
			for x in X:
				if type(x["value"][0]) != G1Elem:
					print("element in C not from G1", type(x["value"][0]))
					return 0
				if "committed" in x:
					if x["committed"][1] != G1Elem.inf(self.G) or x["committed"][2] != G2Elem.inf(self.G):
						print("not a commitment of Zp")
						return 0
			for y in Y:
				if type(y["value"][0]) != G2Elem:
					print("element in D not from G2", type(y["value"][0]))
					return 0
				if "committed" in y:
					if y["committed"][1] != G1Elem.inf(self.G) or y["committed"]["Zp"] != 0:
						print("not a commitment of G2")
						return 0
		elif eq in ["QE", "QC1", "QC2"]:
			for x in X:
				if type(x["value"][0]) != G1Elem :
					print("element in C not from Zp", type(x["value"][0]))
					return 0
				if "committed" in x:
					if x["committed"][1] != G1Elem.inf(self.G) or x["committed"][2] != G2Elem.inf(self.G):
						print("not a commitment of Zp")
						return 0
			for y in Y:
				if type(y["value"][0]) != G2Elem:
					print("element in D not from Zp", type(y["value"][0]))
					return 0
				if "committed" in y:
					if y["committed"][1] != G1Elem.inf(self.G) or y["committed"][2] != G2Elem.inf(self.G):
						print("not a commitment of Zp")
						return 0
		else:
			print("ERROR: unknown eq", eq)
			return 0

		for i in range(len(X)):
			for j in range(len(Y)):
				if type(C[i][j]) != Bn and type(C[i][j]) != int:
					print("element in Gamma not from Zp", i,j, C[i][j], type(C[i][j]))
					return 0
					
		return 1

	def CheckEqResult(self, eq, T):
		#print("CheckEqResult")
		if eq in ["PPE", "PN1", "PC1", "PN2", "PC2"]:
			if T != GTElem.zero(self.G):
				return 0
		elif eq in ["ME1", "MN1", "MC1", "ML1"]:
			if T != G1Elem.inf(self.G):
				return 0
		elif eq in ["ME2", "MN2", "MC2", "ML2"]:
			if T != G2Elem.inf(self.G):
				return 0
		elif eq in ["QE", "QC1", "QC2"]:
			if T != Bn(0):
				return 0

		return 1

	def CheckEqRelation(self, eq, X, Y, C, T):
		#print("CheckEqRelation")
		if self.CheckEqType(eq, X, Y, C) == 0:
			return 0
		if self.CheckEqResult(eq, T) == 0:
			return 0
		return 1


	def MakeProof(self, X, Y, C):
		#print("MakeProof")
		pi2_v1 = [G2Elem.inf(self.G), G2Elem.inf(self.G)]
		pi2_w1 = [G2Elem.inf(self.G), G2Elem.inf(self.G)]

		pi1_v2 = [G1Elem.inf(self.G), G1Elem.inf(self.G)]
		pi1_w2 = [G1Elem.inf(self.G), G1Elem.inf(self.G)]

		#count_exp_g1 = 0
		#count_add_g1 = -1
		#count_exp_g2 = 0
		#count_add_g2 = -1
		
		for i in range(len(X)):
			xi = X[i]["value"]
			xr = 0
			xs = 0
			if "committed" in X[i]:
				xr = X[i]["committed"]["r"]
				xs = X[i]["committed"]["s"] 
			for j in range(len(Y)):
				yj = Y[j]["value"]
				yr = 0
				ys = 0
				if "committed" in Y[j]:
					yr = Y[j]["committed"]["r"]
					ys = Y[j]["committed"]["s"]

				cij = C[i][j]

				if cij != 0:
					for vec in range(2):
						""" We need to work on the commitment of x|y and not on the value of x|y not to make any assumption on the type of x|y"""
						if  xr != 0:
							pi2_v1[vec] += (xr * cij) * yj[vec]
							#count_exp_g1 += 1
							#count_add_g1 += 1
						if xs != 0:
							pi2_w1[vec] += (xs * cij) * yj[vec]
							#count_exp_g1 += 1
							#count_add_g1 += 1
						if yr != 0:
							temp = xi[vec]
							if xr != 0:
								temp = temp - self.v1[vec] * xr
								#count_exp_g2 += 1
								#count_add_g2 += 1
							if xs != 0:
								temp = temp - self.w1[vec] * xs
								#count_exp_g2 += 1
								#count_add_g2 += 1
							pi1_v2[vec] += temp * (cij * yr)
							#count_exp_g2 += 1
							#count_add_g2 += 1
						if ys != 0:
							temp = xi[vec]
							if xr != 0:
								temp = temp - self.v1[vec] * xr
								#count_exp_g2 += 1
								#count_add_g2 += 1
							if xs != 0:
								temp = temp - self.w1[vec] * xs
								#count_exp_g2 += 1
								#count_add_g2 += 1
							pi1_w2[vec] += temp * (cij * ys)
							#count_exp_g2 += 1
							#count_add_g2 += 1
		#print("Exp in G1", count_exp_g1)
		#print("Exp in G2", count_exp_g2)
		#print("Add in G1", count_add_g1)
		#print("Add in G2", count_add_g2)
		return pi2_v1, pi2_w1, pi1_v2, pi1_w2

	def Randomize(self, eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2):
		#print("Randomize")
		a = Bn(0)
		b = Bn(0)
		c = Bn(0)
		d = Bn(0)

		#count_exp_g1 = 0
		#count_add_g1 = -1
		#count_exp_g2 = 0
		#count_add_g2 = -1

		if eq == "PPE":
			a = self.order.random()
			b = self.order.random()
			c = self.order.random()
			d = self.order.random()

		elif eq in ["PN1", "ME2"]: 
			a = self.order.random()
			b = self.order.random()

		elif eq in ["PN2", "ME1"]: 
			a = self.order.random()
			c = self.order.random()

		elif eq in ["MN1", "MN2", "QE"]: 
			a = self.order.random()

		for i in range(2):
			if a != 0:
				pi2_v1[i] += a*self.v2[i]
				pi1_v2[i] += -a*self.v1[i]
				#count_exp_g1 += 1
				#count_exp_g2 += 1
				#count_add_g1 += 1
				#count_add_g2 += 1
			if b != 0:
				pi2_v1[i] += b*self.w2[i]
				pi1_w2[i] += -b*self.v1[i]
				#count_exp_g1 += 1
				#count_exp_g2 += 1
				#count_add_g1 += 1
				#count_add_g2 += 1
			if c != 0:
				pi2_w1[i] += c*self.v2[i]
				pi1_v2[i] += -c*self.w1[i]
				#count_exp_g1 += 1
				#count_exp_g2 += 1
				#count_add_g1 += 1
				#count_add_g2 += 1
			if d != 0:
				pi2_w1[i] += d*self.w2[i]
				pi1_w2[i] += -d*self.w1[i]
				#count_exp_g1 += 1
				#count_exp_g2 += 1
				#count_add_g1 += 1
				#count_add_g2 += 1
			
			#pi2_v1[i] +=  a*self.v2[i] + b*self.w2[i]
			#pi2_w1[i] +=  c*self.v2[i] + d*self.w2[i]
			#pi1_v2[i] += -a*self.v1[i] - c*self.w1[i]
			#pi1_w2[i] += -b*self.v1[i] - d*self.w1[i]

		#print("Exp in G1", count_exp_g1)
		#print("Exp in G2", count_exp_g2)
		#print("Add in G1", count_add_g1)
		#print("Add in G2", count_add_g2)

		return pi2_v1, pi2_w1, pi1_v2, pi1_w2

	def MakeMatrices(self, x, b, a, y, c, t):
		#print("Make Matrices")
		X = []
		X.extend(x)
		X.extend(a)

		Y = [] 
		Y.extend(y)
		Y.extend(b)

		C = []
		for i in range(len(x)+len(a)):
			row = []
			for j in range(len(y)+ len(b)):
				if i < len(x):
					if j < len(y):
						temp = c[i][j]
						if type(temp) != Bn:
							temp = Bn(temp)
						row.append(c[i][j])
					elif j == len(y) + i:
						temp = Bn(1)
						if b[j-len(y)]["committed"][2] == G2Elem.inf(self.G) and b[j-len(y)]["type"] == "pub":
							temp = Bn(0)
						row.append(temp)
					else:
						row.append(Bn(0))
				else:
					if i == len(x) + j:
						temp = Bn(1)
						if a[i-len(x)]["committed"][1] == G1Elem.inf(self.G) and a[i-len(x)]["type"] == "pub":
							temp = Bn(0)
						row.append(temp)
					else:
						row.append(Bn(0))
			C.append(row)

		return X, Y, C, t

	def Prove_eq(self, eq, x, b, a, y, c, t):
		#print("Prove_eq")
		X, Y, C, T = self.MakeMatrices(x, b, a, y, c, t)

		return self.Prove(eq, X, Y, C, T)


	def Prove(self, eq, X, Y, C, T):
		#print("Prove", eq, len(X), len(Y), len(C))

		if self.CheckEqFormat(eq, X, Y, C) == 0:
			print("Proof Error: check Format")
			return 0, []
		if self.CheckEqRelation(eq, X, Y, C, T) == 0:
			print("Proof Error: check Relation")
			return 0, []

		pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = self.MakeProof(X, Y, C)
		pi2_v1, pi2_w1, pi1_v2, pi1_w2 = self.Randomize(eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
		X, Y = self.CleanMatrices(X, Y)

		return 1, [eq, X, Y, C, T, pi2_v1, pi2_w1, pi1_v2, pi1_w2]

	def Prove_aggreg(self, eq,x, b, a, y, c, t):
		#print("Prove", eq, len(X), len(Y), len(C))
		X, Y, C, T = self.MakeMatrices(x, b, a, y, c, t)

		if self.CheckEqFormat(eq, X, Y, C) == 0:
			print("Proof Error: check Format")
			return 0, []
		if self.CheckEqRelation(eq, X, Y, C, T) == 0:
			print("Proof Error: check Relation")
			return 0, []

		pi2_v1, pi2_w1, pi1_v2, pi1_w2 = self.MakeProof(X, Y, C)
		X, Y = self.CleanMatrices(X, Y)

		return 1, [eq, X, Y, C, T, pi2_v1, pi2_w1, pi1_v2, pi1_w2]

	def CleanMatrices(self, X, Y):
		for i in range(len(X)):
			x = {"type":X[i]["type"], "value":X[i]["value"]}
			X[i]=x
		for j in range(len(Y)):
			y = {"type":Y[j]["type"], "value":Y[j]["value"]}
			Y[j]=y
		return X, Y

	def CommitProof_eq(self, eq, x, b, a, y, c, t):
		#print("CommitProof_eq")

		if not self.verifyEq(eq, x, b, a, y, c, t):
			print("uncommitted eq does not verify")
			return 0, []

		try:
			for i in range(len(x)):
				x[i] = self.CommitG1(x[i])
			for i in range(len(a)):
				a[i] = self.CommitG1(a[i])

			for j in range(len(y)):
				y[j] = self.CommitG2(y[j])
			for j in range(len(b)):
				b[j] = self.CommitG2(b[j])

			return self.Prove_eq(eq, x, b, a, y, c, t)
		except Exception as e:
			print("Error:", e)


	def verifyEq(self, eq, x, b, a, y, c, t):
		#print("verifyEq")
		if eq in ["PPE", "PN1", "PC1", "PN2", "PC2"]:
			#print("eq in [\"PPE\", \"PN1\", \"PC1\", \"PN2\", \"PC2\"]")
			T = GTElem.zero(self.G)
			for i in range(min(len(x), len(b))):
				T = T * self.e(x[i]["value"], b[i]["value"])
			for j in range(min(len(a),len(y))):
				T = T * self.e(a[j]["value"], y[j]["value"])
			for i in range(len(c)):
				for j in range(len(c[i])):
					T = T * self.e(c[i][j]*x[i]["value"], y[j]["value"])
			return T.eq(t)
		else :
			#print("eq NOT in [\"PPE\", \"PN1\", \"PC1\", \"PN2\", \"PC2\"]")
			T = Bn(0)
			if eq in ["ME1", "MN1", "MC1", "ML1"]:
				T = G1Elem.inf(self.G)
			elif eq in ["ME2", "MN2", "MC2", "ML2"]:
				T = G2Elem.inf(self.G)
			elif eq not in ["QE", "QC1", "QC2"]:
				print("eq error", eq)
				return 0
			for i in range(min(len(x), len(b))):
				T += x[i]["value"] * b[i]["value"]
			for j in range(min(len(a),len(y))):
				T += a[j]["value"] * y[j]["value"]
			for i in range(len(c)):
				for j in range(len(c[i])):
					if c[i][j] != 0:				
						T += c[i][j] * x[i]["value"] * y[j]["value"]
			return T.eq(t)


	def CommitProof(self, eq, X, Y, C, T):
		#print("CommitProof")
		try:
			for i in range(len(X)):
				X[i] = self.CommitG1(X[i])

			for j in range(len(Y)):
				Y[j] = self.CommitG2(Y[j])

			return self.Prove(eq, X, Y, C, T)
		except Exception as e:
			print("Error:", e)

	def Verify(self, eq, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2):
		#print("Verify")
		if self.CheckEqFormat(eq, X, Y, C) == 0:
			print("Verify: Error - eqformat")
			return 0

		for i in range(2):
			if type(pi2_v1[i]) != G2Elem:
				print("Verify: Error - pi2v1 type")
				return 0
			if type(pi2_w1[i]) != G2Elem:
				print("Verify: Error - pi2w1 type")
				return 0
			if type(pi1_v2[i]) != G1Elem:
				print("Verify: Error - pi1v2 type")
				return 0
			if type(pi1_w2[i]) != G1Elem:
				print("Verify: Error - pi1w2 type")
				return 0
		
		res_eq = [GTElem.one(self.G), GTElem.one(self.G)]

		#count_exp_g1 = 0
		#pairing = -1
		#count_add_gt = -1
		
		for i in range(len(X)):
			for j in range(len(Y)):
				cij = C[i][j]
				xi = X[i]["value"]
				yj = Y[j]["value"]
				if cij != 0 :
					for vec in range(2):
						res_eq[vec] = res_eq[vec] * self.e(cij*xi[vec], yj[vec])
						#count_exp_g1 += 1
						#count_add_gt += 1
						#pairing += 1

		res_proof = [GTElem.one(self.G), GTElem.one(self.G)]
		for vec in range(2):
			res_proof[vec] = res_proof[vec] * self.e(self.v1[vec], pi2_v1[vec])
			res_proof[vec] = res_proof[vec] * self.e(self.w1[vec], pi2_w1[vec])
			res_proof[vec] = res_proof[vec] * self.e(pi1_v2[vec], self.v2[vec])
			res_proof[vec] = res_proof[vec] * self.e(pi1_w2[vec], self.w2[vec])
			#count_add_gt += 4
			#pairing += 4

		#print("Exp G1", count_exp_g1)
		#print("Add GT", count_add_gt)
		#print("Pairing", pairing)

		res = 1
		for vec in range(2):
			if res_eq[vec] != res_proof[vec]:
				print("Verify: Equation", vec, "does not verify")
				res = 0

		return res

	
	def Calc_PPE(self, x1, b2, a1, y2, c, t):
		n = len(x1)
		m = len(y2)
		if len(b2) != n or len(a1) != m or len(c) != n*m or len(t) != 1:
			return 0

		Bt = GTElem.one(self.G)
		for i in range(m):
			Bt = Bt * self.e(x1[i], b2[i])

		At = GTElem.one(self.G)
		for j in range(n):
			At = At*self.e(a1[j],y2[j])

		Ct = GTElem.one(self.G)
		for i in range(n):
			for j in range(m):
				Ct = Ct*self.e(c[i][j]*x1[i],y2[j])

		validate = 0
		if At*Bt*Ct == t:
			validate = 1
		return At*Bt*Ct, t, validate

	def Calc_MSG1(self, x1, b, a1, y, c, t1):
		n = len(x1)
		m = len(y)
		if len(b) != n or len(a1) != m or len(c) != n*m or len(t1) != 1:
			return 0

		B1 = G1Elem.inf(self.G)
		for i in range(m):
			B1 += b[i] * x1[i]

		A1 = G1Elem.inf(self.G)
		for j in range(n):
			A1 += y[j] * a1[j]

		C1 = G1Elem.inf(self.G)
		for i in range(n):
			for j in range(m):
				C1 += (c[i][j] * y[j]) * x1[i]

		validate=0
		if A1+B1+C1 == t1:
			validate = 1
		return A1+B1+C1, t1, validate

	def Calc_MSG2(self, x, b2, a, y2, c, t2):
		n = len(x)
		m = len(y2)
		if len(b2) != n or len(a) != m or len(c) != n*m or len(t2) != 1:
			return 0
		B2 = G2Elem.one(self.G)
		for i in range(m):
			B2 += x[i] * b2[i]

		A2 = G2Elem.one(self.G)
		for j in range(n):
			A += a[j] * y2[j]

		C2=G2Elem.one(self.G)
		for i in range(n):
			for j in range(m):
				C2+= (c[i][j] * x[j]) * y2[i]

		validate = 0
		if A2+B2+C2 == t2:
			z=1
		return A2+B2+C2, t2, validate

	def Calc_QUAD(self, x, b, a, y, c, t):
		n = len(x)
		m = len(y)
		if len(b) != n or len(a) != m or len(c) != n*m or len(t) != 1:
			return 0

		B = Bn(0)
		for i in range(m):
			B += x[i] * b[i]

		A = Bn(0)
		for j in range(n):
			A += a[j] * y[j]

		C=Bn(0)
		for i in range(n):
			for j in range(m):
				C += (c[i][j] * x[j]) * y[i]

		validate = 0
		if A+B+C == t:
			validate=1
		return A+B+C, t , validate



	def proof_usps_hidesigandsigner(
			self,
			M=[b"Hello World!", b"Hello me!", b"Hello us!"]
		):
		"""" creates GS proof that sig verifies with pk, signature and first message secret"""
		#print("SPS: Prove")

		self.ExtGen()
		params = (self.G, self.order, self.g1, self.g2, self.e)

		sps = USPS()
		sk, pk = sps.keygen(params)
		gz, gr, pki = pk

		m = []
		for i in range(len(M)):
			if type(M[i]) == bytes:
				m.append(self.G.hashG1(M[i]))
			elif type(M[i]) == G1Elem:
				m.append(M[i])
			else:
				 return 0, []

		if len(m) < sps.n:
			for i in range(sps.nb_msg - len(m)):
				m.append(G1Elem.inf(self.G))
		elif sps.nb_msg < len(m):
			return
		sig = sps.sign(params, sk,  m)

		if sps.verify(params, pk, m, sig) ==0:
			print("Signature does not verify")
			return 

		print(len(m))
		X = [{"type":"com", "value":sig[0]}, {"type":"com", "value":sig[1]}]
		B = [{"type":"pub", "value":G2Elem.inf(self.G)}, {"type":"pub", "value":G2Elem.inf(self.G)}]

		A = [{"type":"pub", "value":G1Elem.inf(self.G)}, {"type":"pub", "value":G1Elem.inf(self.G)}]
		Y = [{"type":"com", "value":gz}, {"type": "com", "value":gr}]

		for i in range(len(m)):
			if i == 0:
				X.append({"type":"pub", "value":m[i]})
			else:
				X.append({"type":"com", "value":m[i]})
			B.append( {"type":"pub", "value":G2Elem.inf(self.G)} )

		for j in range(len(pki)):
			Y.append({"type": "com", "value":pki[j]})
			A.append( {"type":"pub", "value":G1Elem.inf(self.G)} )

		C = []
		for i in range(len(X)):
			row = []
			for j in range(len(Y)):
				var = Bn(0)
				if i == j:
					var = Bn(1)
				row.append(var)
			C.append(row)
		print(C)

		success, res = self.CommitProof_eq("PPE", X, B, A, Y, C, GTElem.zero(self.G))
		verify = 0
		if success:
			eq_type, X1, Y1, C1, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res
			pi2_v1, pi2_w1, pi1_v2, pi1_w2 = self.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
			verify = self.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
			if verify:
				res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
		print(success, verify)
		print()

		return verify, res
		
	def proof_bsps_hidemandsig(
		self,
		M1 = [b"Hello World!", b"Hello me!", b"Hello us!"], 
		M2 = [b"Hello you!"]):
		""" create GS proof that sig verifies with the signature and all but the first message secret """

		params = (self.G, self.order, self.g1, self.g2, self.e)

		sps = BSPS()
		sk, pk = sps.keygen(params)

		u, w, v, z = sk
		U, W, V, Z = pk

		m1 = []
		for i in range(len(M1)):
			if type(M1[i]) == bytes:
				m1.append(self.G.hashG1(M1[i]))
			elif type(M1[i]) == G1Elem:
				m1.append(M1[i])
			else:
				 return 0, []

		m2 = []
		for j in range(len(M2)):
			if type(M2[j]) == bytes:
				m2.append(Bn.from_binary(M2[j]).mod_add(Bn(0), self.order) * self.g2)
			elif type(M2[j]) == G2Elem:
				m2.append(M2[j])
			else:
				return 0, []

		if len(m1) < sps.msg_g1:
			for i in range(sps.msg_g1 - len(m1)):
				m1.append(G1Elem.inf(self.G))
		elif sps.msg_g1 < len(m1):
			return 0, []
		
		if len(m2) < sps.msg_g2:
			for i in range(sps.msg_g2 - len(m2)):
				m2.append(G2Elem.inf(self.G))
		elif sps.msg_g2 < len(m2):
			return 0, []

		res = []

		R, S, T = sps.sign(params, sk, m1, m2)
		print("signature verifies?", sps.verify(params, pk, m1, m2, (R,S,T)))

		print("first equation")
		x1 = [{"type":"com", "value":R}, {"type":"com", "value":S}, {"type":"bas", "value":self.g1}]
		b1 = [{"type":"pub", "value":V}, {"type":"bas", "value":self.g2}, {"type":"pub", "value":Z*Bn(-1)}]
		a1 = []
		y1 = []

		#res1 = e(R,V) * e(S,g2)
		#res1.eq(e(g1,Z))
		for i in range(len(m1)):
			if i == 0:
				x1.append({"type":"pub", "value":m1[i]})
				b1.append({"type":"bas", "value":W[i]})
			else:				
				x1.append({"type":"com", "value":m1[i]})
				b1.append({"type":"pub", "value":W[i]})
			#res1 = res1 * e(m1[i],W[i])

		c1 = []
		for i in range(len(x1)):
			row = []
			for j in range(len(y1)):
				c = Bn(0)
				row.append(c)
			c1.append(row)
		print(c1)


		success1, res1 = self.CommitProof_eq("PPE", x1, b1, a1, y1, c1, GTElem.zero(self.G))
		verify1 = 0
		if success1:
			eq_type, X1, Y1, C1, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res1
			pi2_v1, pi2_w1, pi1_v2, pi1_w2 = self.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
			verify1 = self.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
			if verify1:
				res1 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
				res.append(res1)
		print(success1, verify1)
		print()

		print("second equation")
		x2 = [{"type":"com", "value":R}, 					{"type":"bas", "value":self.g1}]
		b2 = [{"type":"pub", "value":G2Elem.inf(self.G)},	{"type":"pub", "value":G2Elem.inf(self.G)}]
		a2 = [{"type":"pub", "value":G1Elem.inf(self.G)},	{"type":"pub", "value":G1Elem.inf(self.G)}]
		y2 = [{"type":"com", "value":T}, 					{"type":"bas", "value":self.g2}]

		#res2 = e(R,T)
		#res2.eq(e(g1,g2))
		for j in range(len(m2)):
			a2.append({"type":"pub", "value":U[j]})
			y2.append({"type":"com", "value":m2[j]})
			#res2 = res2 * e(U[j],m2[j])

		c2 = []
		for i in range(len(x2)):
			row = []
			for j in range(len(y2)):
				c = Bn(0)
				if (i == 0 and j == 0):
					c = Bn(1)
				if (i == 1 and j == 1):
					c = Bn(-1)
				row.append(c)
			c2.append(row)
		print(c2)

		
		success2, res2 = self.CommitProof_eq("PPE", x2, b2, a2, y2, c2, GTElem.zero(self.G))
		verify2 = 0
		if success2:
			eq_type, X2, Y2, C2, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res2
			pi2_v1, pi2_w1, pi1_v2, pi1_w2 = self.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
			verify2 = self.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
			if verify2:
				res2 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X2, Y2, C2, T_eq] ]
				res.append(res2)
		print(success2, verify2)
		print()

		print("success?", success1, success2)
		print("verify?", verify1, verify2)

		return verify1*verify2, res

def proof_exponent_gs1_all_secret(gsp, g1, sk, pk):
	x = [{"type":"com", "value":g1}, {"type":"com", "value":pk}]
	b = [{"type":"sca", "value":Bn(0)}, {"type":"unt", "value":Bn(1)}]
	a = [{"type":"pub", "value":G1Elem.inf(gsp.G)}]
	y = [{"type":"sca", "value":sk}]
	c = [[Bn(-1)],[Bn(0)]]

	success, res = gsp.CommitProof_eq("ME1", x, b, a, y, c, G1Elem.inf(gsp.G))
	verify = 0
	if success:
		eq_type, X, Y, C, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res
		pi2_v1, pi2_w1, pi1_v2, pi1_w2 = gsp.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
		verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
		if verify:
			res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X, Y, C, T_eq] ]
	print("Do we successfully create a first proof?", success)
	print("Does the first proof successfully verify?", verify)
	print()
	return verify, res

def proof_exponent_gs1_pk_public(gsp, g1, sk, pk):
	x = [{"type":"com", "value":g1}, {"type":"com", "value":pk}]
	b = [{"type":"sca", "value":Bn(0)}, {"type":"unt", "value":Bn(1)}]
	a = [{"type":"pub", "value":G1Elem.inf(gsp.G)}]
	y = [{"type":"sca", "value":sk}]
	c = [[Bn(-1)],[Bn(0)]]

	success, res = gsp.CommitProof_eq("ME1", x, b, a, y, c, G1Elem.inf(gsp.G))
	verify = 0
	if success:
		eq_type, X, Y, C, T_eq, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap = res
		pi2_v1, pi2_w1, pi1_v2, pi1_w2 = gsp.Randomize(eq_type, pi2_v1_ap, pi2_w1_ap, pi1_v2_ap, pi1_w2_ap)
		verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
		if verify:
			res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X, Y, C, T_eq] ]
	print("Do we successfully create a first proof?", success)
	print("Does the first proof successfully verify?", verify)
	print()
	return verify, res

		
# ---------- TESTS -------------

def test_ExtGen():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()


def test_G1_group_comm():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck
	x1 = order.random()*G.gen1()
	r = order.random()
	s = order.random()
	struct = {}
	struct = GS.CommitG1Group(struct, x1, r, s)

	assert struct["value"][0] == G1Elem.inf(G) + r * GS.v1[0] + s * GS.w1[0]
	assert struct["value"][1] == x1 + r * GS.v1[1] + s * GS.w1[1]

	assert struct["committed"][1] == x1
	assert struct["committed"][2] == G2Elem.inf(G)
	assert struct["committed"]["Zp"] == Bn(0)
	assert struct["committed"]["r"] == r
	assert struct["committed"]["s"] == s


def test_G1_scalar_comm():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck
	x = order.random()
	r = order.random()
	struct = {}
	struct = GS.CommitG1Scalar(struct, x, r)

	assert struct["value"][0] == x * GS.u1[0] + r * GS.v1[0]
	assert struct["value"][1] == x * GS.u1[1] + r * GS.v1[1]

	assert struct["committed"][1] == G1Elem.inf(G)
	assert struct["committed"][2] == G2Elem.inf(G)
	assert struct["committed"]["Zp"] == x
	assert struct["committed"]["r"] == r
	assert struct["committed"]["s"] == Bn(0)


def test_G2_group_comm():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck
	y2 = order.random()*G.gen2()
	r = order.random()
	s = order.random()
	struct = {}
	struct = GS.CommitG2Group(struct, y2, r, s)

	assert struct["value"][0] == G2Elem.inf(G) + r * GS.v2[0] + s * GS.w2[0]
	assert struct["value"][1] == y2 + r * GS.v2[1] + s * GS.w2[1]

	assert struct["committed"][1] == G1Elem.inf(G)
	assert struct["committed"][2] == y2
	assert struct["committed"]["Zp"] == Bn(0)
	assert struct["committed"]["r"] == r
	assert struct["committed"]["s"] == s

def test_G2_scalar_comm():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck
	y = order.random()
	r = order.random()
	struct = {}
	struct = GS.CommitG2Scalar(struct, y, r)

	assert struct["value"][0] == y * GS.u2[0] + r * GS.v2[0]
	assert struct["value"][1] == y * GS.u2[1] + r * GS.v2[1]

	assert struct["committed"][1] == G1Elem.inf(G)
	assert struct["committed"][2] == G2Elem.inf(G)
	assert struct["committed"]["Zp"] == y
	assert struct["committed"]["r"] == r
	assert struct["committed"]["s"] == Bn(0)

def test_exp_g1_all_secret():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, g1, g2, e = P
	h = order.random() * g1
	sk = order.random()
	pk = sk * h

	verify, res = proof_exponent_gs1_all_secret(GS, h, sk, pk)
	assert verify

def test_exp_g1_pk_public():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, g1, g2, e = P
	h = order.random() * g1
	sk = order.random()
	pk = sk * h

	verify, res = proof_exponent_gs1_pk_public(GS, h, sk, pk)
	assert verify

def test_ppe_wrongformat():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck

	x = []
	y = []
	a = []
	b = [] 
	for xtype in GS.eqs["PPE"][1]:
		for ytype in GS.eqs["PPE"][2]:
			x1 = order.random() * GS.g1
			y2 = order.random() * GS.g2
			x.append( {"type": "pub", "value": x1} )
			a.append( {"type": ytype, "value": x1} )
			y.append( {"type": ytype, "value": y2} )
			b.append( {"type": "pub", "value": y2} )

	c = []
	for i in range(len(x)):
		row = []
		for j in range(len(y)):
			if i == j :
				row.append(Bn(-2))
			else :
				row.append(0)
		c.append(row)

	t = GTElem.zero(G)

	try :
		success, res = GS.CommitProof_eq("PPE", x, b, a, y, c, t)
	except Exception as e:
		print(e)
	else :
		assert success == 0


def test_all():
	GS = GSProof()
	ck, xk, P = GS.ExtGen()
	order, G, u1, v1, w1, u2, v2, w2 = ck

	for eq_type in GS.eqs.keys():
		print("=============================== ",eq_type)
		x = []
		y = []
		a = []
		b = [] 
		for xtype in GS.eqs[eq_type][1]:
			for ytype in GS.eqs[eq_type][2]:
				if eq_type != "PPE" or (eq_type == "PPE" and xtype not in ["pub", "enc"] and ytype not in ["pub", "enc"]):
					rand1 = order.random()
					x1 = rand1 * GS.g1
					if xtype in ["bas"]:
						rand1 = 1
						x1 = GS.g1
					if xtype in ["sca"]:
						x1 = rand1
					elif xtype in ["unt"]:
						x1 = 1
					rand2 = order.random()
					y2 = rand2 * GS.g2
					if ytype in ["bas"]:
						rand2 = 1
						y2 = GS.g2
					if ytype in ["sca"]:
						y2 = rand2
					elif ytype in ["unt"]:
						y2 = 1
					#print(xtype, rand1, x1)
					#print(ytype, rand2, y2)
					x.append( {"type": xtype, "value": x1} )
					a.append( {"type": xtype, "value": x1} )
					y.append( {"type": ytype, "value": y2} )
					b.append( {"type": ytype, "value": y2} )

		c = []
		for i in range(len(x)):
			row = []
			for j in range(len(y)):
				if i == j :
					row.append(Bn(-2))
				else :
					row.append(Bn(0))
			c.append(row)

		t = GTElem.zero(G)
		if eq_type in ["ME1", "MN1", "MC1", "ML1"]:
			t = G1Elem.inf(G)
		if eq_type in ["ME2", "MN2", "MC2", "ML2"]:
			t = G2Elem.inf(G)
		if eq_type in ["QE", "QC1", "QC2"]:
			t = 0

		try :
			success, res = GS.CommitProof_eq(eq_type, x, b, a, y, c, t)
			print("success ?", success)
		except Exception as e:
			print(e)
		else :
			if success :
				eq, X, Y, C, T, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res
				print("verify ?", GS.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2))
				assert(GS.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2))
	
