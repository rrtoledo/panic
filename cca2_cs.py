from bplib.bp import BpGroup
from bplib.bp import G1Elem
from hashlib import sha256
from petlib.bn import Bn
import pytest

class CSEnc:

    def __init__(self):
        print("CSEnc: Init")

    def params_gen(self):
        """Generates the AHEG for an EC group nid"""
        print("CCA2EG: Setup")
        G = BpGroup()
        g1 = G.gen1()
        o =  G.order()
        g2 = o.random() * g1
        return [o, G, g1, g2]

    def key_gen(self, params):
        """Generates a fresh key pair"""
        print("CSEnc: Key Gen")
        o, _, g1, g2 = params
        x1 = o.random()
        x2 = o.random()
        y1 = o.random()
        y2 = o.random()
        z = o.random()
        priv = [x1, x2, y1, y2, z]
        c = x1 * g1 + x2 * g2
        d = y1 * g1 + y2 * g2
        h =  z * g1
        pub = [c, d, h]
        return [pub, priv]

    def enc(self, params, pub, m):
        """Encrypts the values of a group element"""
        #print("CSEnc: Enc")
        o, G, g1, g2 = params
        c, d, h = pub

        r = o.random()
        u1, u2 = r*g1, r*g2
        e = r*h+m
        
        H = sha256()
        H.update(u1.export() + u2.export() + e.export())
        a = Bn.from_hex(H.hexdigest())

        v = r*c + (a*r)*d

        return [u1, u2, e, v], [r, a]


    def dec(self, priv, cipher):
        """Decrypts an encrypted group element"""
        #print("CSEnc: Dec")
        u1, u2, e, v = cipher
        x1, x2, y1, y2, z = priv

        H = sha256()
        H.update(u1.export() + u2.export() + e.export())
        a = Bn.from_hex(H.hexdigest())

        if v.eq((x1+y1*a)*u1 + (x2+y2*a)*u2):
            return 1, e + (-1*z)*u1
        else:
            print("Decryption error")
            return 0, -1


def test_enc():
    enc = CSEnc()
    params = enc.params_gen()
    o, _, g1, g2 = params
    (pub, priv) = enc.key_gen(params)

    # Check encryption and decryption
    m = o.random() * g1
    c1, _ = enc.enc(params, pub, m)
    b, m1 = enc.dec(priv, c1) 
    print (b)
    assert b == 1 and m1.eq(m)

