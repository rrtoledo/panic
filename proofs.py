from bplib.bp import G1Elem
from bplib.bp import G2Elem
from bplib.bp import GTElem
from petlib.bn import Bn
from gsproof import GSProof
from bsps import BSPS
from usps import USPS
from cca2_cs import CSEnc

def proof_usps_hidesigandsigner(gsp, pk, M, sig):
    """" creates GS proof that a USPS signature verifies
    with the verifying key, the signature and the first message secret"""
    params = gsp.P

    gz, gr = pk[0], pk[1]
    pki = pk[2:]

    X = [sig[0], sig[1]]
    B = []
    A = []
    Y = [gz, gr]

    for i in range(len(M)):
        X.append(M[i])
        Y.append(pki[i])

    C = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            var = Bn(0)
            if i == j:
                var = Bn(1)
            row.append(var)
        C.append(row)

    success, res = gsp.Prove_eq("PPE", X, B, A, Y, C, GTElem.zero(gsp.G))
    verify = 0
    if success:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res
        verify = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify:
            res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
    print("Do we successfully create a proof?", success)
    print("Does the proof successfully verify?", verify)
    return verify, res


def proof_bsps_hidemandsig(gsp, pk, M1, M2, t, sig):
    """ create GS proof that sig verifies with the signature and all but the first message secret """

    params = gsp.P

    U, W, V, Z = pk
    R, S, T = sig

    res = []

    #print("----- first equation")
    x1 = [gsp.Commit({"group":1,"type":"bas", "value":gsp.g1}), R, S]
    b1 = []
    a1 = []
    y1 = [Z, V, gsp.Commit({"group":2, "type":"bas", "value":gsp.g2})]

    #res1 = e(R,V) * e(S,g2)
    #res1.eq(e(g1,Z))
    for i in range(len(W)):
        x1.append(M1[i])
        y1.append(W[i])
        #res1 = res1 * e(m1[i],W[i])

    c1 = []
    for i in range(len(x1)):
        row = []
        for j in range(len(y1)):
            c = Bn(0)
            if i==j:
                c = Bn(1)
                if i == 0:
                    c = Bn(-1)
                if i == 3:
                    if type(t) == int:
                        t = Bn(t)
                    c = t      
            row.append(c)
        c1.append(row)

    success1, res1 = gsp.Prove_eq("PPE", x1, b1, a1, y1, c1, GTElem.zero(gsp.G))
    verify1 = 0
    if success1:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res1
        verify1 = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify1:
            res1 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
            res.append(res1)
    #print("Do we successfully create a first proof?", success1)
    #print("Does the first proof successfully verify?", verify1)
    #print("----- second equation")
    x2 = [gsp.Commit({"group":1,"type":"bas", "value":gsp.g1}), R]
    b2 = []
    a2 = []
    y2 = [gsp.Commit({"group":2, "type":"bas", "value":gsp.g2}), T]

    #res2 = e(R,T)
    #res2.eq(e(g1,g2))
    for j in range(len(U)):
        x2.append(U[j])
        y2.append(M2[j])
        #res2 = res2 * e(U[j],m2[j])

    c2 = []
    for i in range(len(x2)):
        row = []
        for j in range(len(y2)):
            c = Bn(0)
            if i == j:
                c = Bn(1)
                if (i == 0):
                    c = Bn(-1)
            row.append(c)
        c2.append(row)
    
    success2, res2 = gsp.Prove_eq("PPE", x2, b2, a2, y2, c2, GTElem.zero(gsp.G))
    verify2 = 0
    if success2:
        eq_type, X2, Y2, C2, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res2
        verify2 = gsp.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify2:
            res2 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X2, Y2, C2, T_eq] ]
            res.append(res2)
    #print("Do we successfully create a second proof?", success2)
    #print("Does the second proof successfully verify?", verify2)
    print("Are all the proofs successfull created?", success1*success2)
    print("Do all the proofs verify?", verify1*verify2)
    return verify1*verify2, res

def enc_proof(gsp, params, pub, da, m, cipher, rand):

    g1enc, g2enc = params
    c, d, h = pub
    u1, u2, e, v = cipher
    r, a = rand

    res = []
    #print("--- first equation")
    x1 = [e, m, h]
    a1 = []
    y1 = [gsp.Commit({"group":2, "type":"unt", "value":1}), r]
    b1 = []
    c1 = [[-1,0],[1,0],[0,1]]

    success1, res1 = gsp.Prove_eq("ME1", x1, b1, a1, y1, c1, G1Elem.inf(gsp.G))
    verify1 = 0
    if success1:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res1
        verify1 = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify1:
            res1 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X1, Y1, C1, T_eq] ]
            res.append(res1)
    #print("Do we successfully create a first proof?", success1)
    #print("Does the first proof successfully verify?", verify1)

    #print("--- second equation")
    x2 = [g1enc, u1]
    b2 = []
    a2 = []
    y2 = [r, gsp.Commit({"group":2, "type":"unt", "value":1})]
    c2 = [[-1,0],[0,1]]
    success2, res2 = gsp.Prove_eq("MC1", x2, b2, a2, y2, c2, G1Elem.inf(gsp.G))
    verify2 = 0
    if success2:
        eq_type, X2, Y2, C2, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res2
        verify2 = gsp.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify2:
            res2 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X2, Y2, C2, T_eq] ]
            res.append(res2)
    #print("Do we successfully create a second proof?", success2)
    #print("Does the second proof successfully verify?", verify2)

    #print("--- third equation")
    x3 = [g2enc, u2]
    b3 = []
    a3 = []
    y3 = [r, gsp.Commit({"group":2, "type":"unt", "value":1})]
    c3 = [[-1,0],[0,1]]
    success3, res3 = gsp.Prove_eq("MC1", x3, b3, a3, y3, c3, G1Elem.inf(gsp.G))
    verify3 = 0
    if success3:
        eq_type, X3, Y3, C3, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res3
        verify3 = gsp.Verify(eq_type, X3, Y3, C3, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify3:
            res3 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X3, Y3, C3, T_eq] ]
            res.append(res3)
    #print("Do we successfully create a third proof?", success3)
    #print("Does the third proof successfully verify?", verify3)

    #print("--- intermediary equation")
    x_int = [da, d]
    b_int = []
    a_int = []
    y_int = [gsp.Commit({"group":2, "type":"unt", "value":1})]
    c_int = [[-1],[a]]
    success_int, res_int = gsp.Prove_eq("MC1", x_int, b_int, a_int, y_int, c_int, G1Elem.inf(gsp.G))
    verify_int = 0
    if success_int:
        eq_type, X_int, Y_int, C_int, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res_int
        verify_int = gsp.Verify(eq_type, X_int, Y_int, C_int, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify_int:
            res_int = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X_int, Y_int, C_int, T_eq] ]
            res.append(res_int)
    #print("Do we successfully create an int. proof?", success_int)
    #print("Does the int. proof successfully verify?", verify_int)

    #print("--- fourth equation")
    x4 = [c, da, v]
    b4 = []
    a4 = []
    y4 = [r, gsp.Commit({"group":2, "type":"unt", "value":1})]
    c4 = [[1,0],[1,0],[0,-1]]
    success4, res4 = gsp.Prove_eq("MC1", x4, b4, a4, y4, c4, G1Elem.inf(gsp.G))
    verify4 = 0
    if success4:
        eq_type, X4, Y4, C4, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res4
        verify4 = gsp.Verify(eq_type, X4, Y4, C4, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify4:
            res4 = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X4, Y4, C4, T_eq] ]
            res.append(res4)
    #print("Do we successfully create a fourth proof?", success4)
    #print("Does the fourth proof successfully verify?", verify4)

    print("Do we successfully create all the proofs?", success1*success2*success3*success4*success_int)
    print("Do all the proofs successfully verify?", verify1*verify2*verify3*verify_int*verify4)
    return verify1*verify2*verify3*verify_int*verify4, res

def proof_exponent(gsp, pk_i, pk_ui, sku):
    x = [pk_i, pk_ui]
    b = []
    a = []
    y = [sku, gsp.Commit({"group":2, "type":"unt", "value":1})]
    c = [[1,0],[0,-1]]

    eq = "MC1"
    if pk_ui["type"] == "com":
        eq = "ME1"
    success, res = gsp.Prove_eq(eq, x, b, a, y, c, G1Elem.inf(gsp.G))
    verify = 0
    if success:
        eq_type, X, Y, C, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = res
        verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        if verify:
            res = [[pi2_v1, pi2_w1, pi1_v2, pi1_w2], [eq_type, X, Y, C, T_eq] ]
    print("Do we successfully create a proof?", success)
    print("Does the proof successfully verify?", verify)
    print()
    return verify, res

def prepare_proofs(auth, vkI, pki, pkv, m, sig_t, t, pk_ui, pk_uv, sku, cipher, ek, r):
    print("Prepare proof")
    print("Prepare proof: Commit")

    import time

    tc = time.time()
    U, W, V, Z = vkI
    cm_U= []
    for i in range(len(U)):
        cm_U.append(auth.GS.Commit({"group":1, "type":"pub", "value":U[i]}))
    cm_W = []
    for i in range(len(W)):
        cm_W.append(auth.GS.Commit({"group":2, "type":"pub", "value":W[i]}))
    cm_V = auth.GS.Commit({"group":2, "type":"pub", "value":V})
    cm_Z = auth.GS.Commit({"group":2, "type":"pub", "value":Z})
    cm_vkI = [cm_U, cm_W, cm_V, cm_Z]

    h_i, vk_i, sig0_i = pki
    cm_h_i = auth.GS.Commit({"group":1, "type":"com", "value":h_i})
    cm_vk_i = []
    for i in range(len(vk_i)):
        cm_vk_i.append(auth.GS.Commit({"group":2, "type":"com", "value":vk_i[i]}))
    cm_sig0_i = []
    for i in range(len(sig0_i)):
        cm_sig0_i.append(auth.GS.Commit({"group":max(i,1), "type":"com", "value":sig0_i[i]}))

    h_v, vk_v, sig0_v = pkv
    cm_h_v = auth.GS.Commit({"group":1, "type":"pub", "value":h_v})

    cm_m = []
    for i in range(len(m)):
        cm_m.append(auth.GS.Commit({"group":1, "type":"pub", "value":m[i]}))

    sig_i, c_it = sig_t
    cm_sig_i = []
    for i in range(len(sig_i)):
        cm_sig_i.append(auth.GS.Commit({"group":1, "type":"com", "value":sig_i[i]}))
    _, sigt_i = c_it
    cm_sigt_i = []
    for i in range(len(sigt_i)):
        cm_sigt_i.append(auth.GS.Commit({"group":max(i,1), "type":"com", "value":sigt_i[i]}))

    #t: used as public scalar constraint (gamma_ij)

    cm_pk_ui = auth.GS.Commit({"group":1, "type":"com", "value":pk_ui})

    cm_pk_uv = auth.GS.Commit({"group":1, "type":"pub", "value":pk_uv})

    cm_sku = auth.GS.Commit({"group":2, "type":"sca", "value":sku})

    cm_cipher = []
    for i in range(len(cipher)):
        cm_cipher.append(auth.GS.Commit({"group":1, "type":"pub", "value":cipher[i]}))
    cm_r = r
    cm_r[0] = auth.GS.Commit({"group":2, "type":"sca", "value":r[0]})

    cm_ek = []
    for i in range(len(ek)):
        cm_ek.append(auth.GS.Commit({"group":1, "type":"pub", "value":ek[i]}))
    cm_da = auth.GS.Commit({"group":1, "type":"pub", "value":auth.ek[1]*r[1]})

    cm_params_enc = []
    cm_params_enc.append(auth.GS.Commit({"group":1, "type":"pub", "value":auth.GS.v1[0]}))
    cm_params_enc.append(auth.GS.Commit({"group":1, "type":"pub", "value":auth.GS.v1[1]}))

    cm_g1 = auth.GS.Commit({"group":1, "type":"bas", "value":auth.GS.g1})
    cm_g2 = auth.GS.Commit({"group":2, "type":"bas", "value":auth.GS.g2})
    cm_1 = auth.GS.Commit({"group":2, "type":"unt", "value":1})

    tcnd = time.time()
    print("--- Commitment time:", tcnd-tc)
    
    print("Prepare proof: Prove")
    tp = time.time()
    print("--- enc proof")
    verify_enc, res_enc = enc_proof(auth.GS, cm_params_enc, cm_ek, cm_da, cm_h_i, cm_cipher, cm_r)

    print("--- sigi proof")
    cm_msg = [cm_pk_ui]
    cm_msg.extend(cm_m)
    for i in range(len(cm_vk_i) - 2 - len(cm_msg)):
        cm_msg.append(auth.GS.Commit({"group":1, "type":"pub", "value":G1Elem.inf(auth.GS.G)}))
    verify_sigi, res_sigi = proof_usps_hidesigandsigner(auth.GS, cm_vk_i, cm_msg, cm_sig_i)

    print("--- sigi0 proof")
    cm_msg1 = [cm_g1, cm_h_i]
    cm_msg2 = cm_vk_i
    verify_sigi0, res_sigi0 = proof_bsps_hidemandsig(auth.GS, cm_vkI, cm_msg1, cm_msg2, Bn(0), cm_sig0_i)

    print("--- sigit proof")
    verify_sigit, res_sigit = proof_bsps_hidemandsig(auth.GS, cm_vkI, cm_msg1, cm_msg2, t, cm_sigt_i)

    print("--- exponent proofs")
    verify_pkui, res_pkui = proof_exponent(auth.GS, cm_h_i, cm_pk_ui, cm_sku)
    verify_pkuv, res_pkuv = proof_exponent(auth.GS, cm_h_v, cm_pk_uv, cm_sku)

    tpend = time.time()
    print("--- Proving time:", tp - tpend)

    verify = verify_enc*verify_sigi*verify_sigi0*verify_sigit*verify_pkui*verify_pkuv
    res = []
    res.extend(res_enc)
    res.append(res_pkui)
    res.append(res_pkuv)
    res.append(res_sigi)
    res.extend(res_sigi0)
    res.extend(res_sigit)
    print("Do all the proofs verify?", verify_enc*verify_sigi*verify_sigi0*verify_sigit*verify_pkui*verify_pkuv)

    return verify, res

def test_enc_proof():
    gsp = GSProof()
    gsp.ExtGen()
    params = gsp.order, gsp.G, gsp.v1[0], gsp.v1[1]

    enc = CSEnc()
    ek, dk = enc.key_gen(params)
    m = gsp.order.random() * gsp.v1[0]
    cipher, rand = enc.enc(params, ek, m)

    Params = []
    Params.append(gsp.Commit({"group":1, "type":"pub", "value":gsp.v1[0]}))
    Params.append(gsp.Commit({"group":1, "type":"pub", "value":gsp.v1[1]}))

    EK = []
    for i in range(len(ek)):
        EK.append(gsp.Commit({"group":1, "type":"pub", "value":ek[i]}))
    da = gsp.Commit({"group":1, "type":"pub", "value":ek[1]*rand[1]})

    M = gsp.Commit({"group":1, "type":"com", "value":m})

    Cipher = []
    for i in range(len(cipher)):
        Cipher.append(gsp.Commit({"group":1, "type":"pub", "value":cipher[i]}))
    
    Rand = [gsp.Commit({"group":2, "type":"sca", "value":rand[0]}), rand[1]]

    verify, res = enc_proof(gsp, Params, EK, da, M, Cipher, Rand)
    assert verify

def test_usps_proof():
    gsp = GSProof()
    gsp.ExtGen()
    params = gsp.P

    sps = USPS()
    sk, pk = sps.keygen(params)
    gz, gr = pk[0], pk[1]
    pki = pk[2:]

    m = [
        gsp.G.hashG1(b"Hello World!"),
        gsp.G.hashG1(b"Hello me!"),
        gsp.G.hashG1(b"Hello us!")]

    if len(m) < len(pk)-2:
        for i in range(len(pk)-2 - len(m)):
            m.append(G1Elem.inf(gsp.G))
    elif len(pk)-2 < len(m):
        return
    sig = sps.sign(params, sk, m)
    print("Does the signature verify?", sps.verify(params, pk, m, sig))

    PK = []
    for i in range(len(pk)):
        PK.append(gsp.Commit({"group":2, "type":"com", "value":pk[i]}))
    M = [
        gsp.Commit({"group":1, "type":"com", "value":m[0]}),
        gsp.Commit({"group":1, "group":1, "type":"pub", "value":m[1]}),
        gsp.Commit({"group":1, "type":"pub", "value":m[2]})]
    for i in range(len(pk)-2 - len(m)):
        M.append(gsp.Commit({"group":1, "type":"pub", "value":G1Elem.inf(gsp.G)}))
    SIG = []
    for i in range(len(sig)):
        SIG.append(gsp.Commit({"group":1, "type":"com", "value":sig[i]}))

    verify, res = proof_usps_hidesigandsigner(gsp, PK, M, SIG)
    assert verify


def test_bsps_proof():
    """ create GS proof that sig verifies with the signature and all but the first message secret """
    gsp = GSProof()
    gsp.ExtGen()
    params = gsp.P

    sps = BSPS()
    sk, pk = sps.keygen(params)

    u, w, v, z = sk
    U, W, V, Z = pk

    t= Bn(7)
    m1 = [gsp.g1*t, gsp.G.hashG1(b"Hello World!")]
    m2 = [Bn.from_binary(b"Hello you!") * gsp.g2]

    for i in range(len(W) - len(m1)):
        m1.append(G1Elem.inf(gsp.G))
    
    for i in range(len(U) - len(m2)):
        m2.append(G2Elem.inf(gsp.G))

    sig = sps.sign(params, sk, m1, m2)
    print("Does the signature verify?", sps.verify(params, pk, m1, m2, sig))

    PK = [[],[]]
    for i in range(len(U)):
        PK[0].append(gsp.Commit({"group":1, "type":"pub", "value":U[i]}))
    for i in range(len(W)):
        PK[1].append(gsp.Commit({"group":2, "type":"pub", "value":W[i]})) 
    PK.append(gsp.Commit({"group":2, "type":"pub", "value":V}))
    PK.append(gsp.Commit({"group":2, "type":"pub", "value":Z}))
    M1 = [gsp.Commit({"group":1,"type":"bas", "value":gsp.g1})]
    for i in range(len(m1)-1):
        M1.append(gsp.Commit({"group":1, "type":"com", "value":m1[i+1]}))
    M2 = []
    for i in range(len(m2)):
        M2.append(gsp.Commit({"group":2, "type":"com", "value":m2[i]}))
    SIG = [
        gsp.Commit({"group":1, "type":"com", "value":sig[0]}),
        gsp.Commit({"group":1, "type":"com", "value":sig[1]}),
        gsp.Commit({"group":2, "type":"com", "value":sig[2]})
    ] 
    verify, result = proof_bsps_hidemandsig(gsp, PK, M1, M2, t, SIG)
    assert verify
