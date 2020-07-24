## An example of the simple Schnorr sigma protocol
## to prove that one knows x, such that h = g^x for 
## a public generator h and g.

from petlib.bn import Bn
from petlib.ec import EcGroup, EcPt

from hashlib import sha256

class SOK:

    def __init__(self):
       print("SOK: Init")

    def challenge(self, elements):
        """Packages a challenge in a bijective way"""
        elem = [len(elements)] + elements
        elem_str = map(str, elem)
        elem_len = map(lambda x: "%s||%s" % (len(x) , x), elem_str)
        state = "|".join(elem_len)
        H = sha256()
        H.update(state.encode("utf8"))
        return H.digest()
        

    def setup(self):
        G = EcGroup(713)
        g = G.generator()
        o = G.order()
        return o, G, g

    def prove(self, params, h, g, x, m=""):
        #print("SOK: Prove")
        """Schnorr proof of the statement ZK(x ; h = g^x)"""
        assert x * g == h
        o, G, _ = params
        w = o.random()
        W = w * g

        state = ['schnorr', o, g.export(), h.export(), m, W.export()]
        hash_c = self.challenge(state)
        c = Bn.from_binary(hash_c) % o
        r = (w - c * x) % o
        return (c, r)

    def verify(self, params, h, g, proof, m=""):
        #print("SOK: Verify")
        """Verify the statement ZK(x ; h = g^x)"""
        o, G, _ = params
        c, r = proof
        W = (r * g + c * h)

        state = ['schnorr', o, g.export(), h.export(), m, W.export()]
        hash_c = self.challenge(state)
        c2 = Bn.from_binary(hash_c) % o
        return c == c2


def test_zkp():
    sok = SOK()
    params = sok.setup()
    o, G, g = params
    x = o.random()
    h = x * g

    ## Use it as a Zk proof
    proof = sok.prove(params, h, g, x)
    assert sok.verify(params, h, g, proof)
    assert not sok.verify(params, g, h, proof)

    ## Use it as a signature scheme
    proofm = sok.prove(params, h, g, x, m = "Hello World!")
    assert sok.verify(params, h, g, proofm, m = "Hello World!")
    assert not sok.verify(params, h, g, proofm, m = "Other String")
