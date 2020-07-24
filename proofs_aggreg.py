from bplib.bp import G1Elem
from bplib.bp import G2Elem
from bplib.bp import GTElem
from petlib.bn import Bn
from gsproof import GSProof
from bsps import BSPS
from usps import USPS
from cca2_cs import CSEnc
from hashlib import sha256

def proof_sigi(gsp, X, Y, M):
    """" creates GS proof that a USPS signature verifies
    with the verifying key, the signature and the first message secret"""
    res = []
    A = []
    B = []
    C = []
    counter_i = len(X)-len(M)-2
    counter_j = len(Y)-len(M)-2
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            var = Bn(0)
            if i == counter_i and j == counter_j:
                var = Bn(1)
                counter_i += 1
                counter_j += 1
            row.append(var)
        C.append(row)
 
    success, result = gsp.Prove_aggreg("PPE", X, B, A, Y, C, GTElem.zero(gsp.G))
    verify = 0
    if success:
        eq_type, X, Y, C, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result
        #verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify:
        res = [C, [pi2_v1, pi2_w1, pi1_v2, pi1_w2] ]
    #print("Do we successfully create a proof?", success)
    #print("Does the proof successfully verify?", verify)

    b = challenge(result)
    for i in range(len(res[0])):
        for j in range(len(res[0][0])):
            if res[0][i][j] != 0:
                res[0][i][j]= b*res[0][i][j]
    for i in range(len(res[1])):
        for j in range(len(res[1][i])):
            res[1][i][j] = b*res[1][i][j]

    #verify = gsp.Verify("PPE", X, Y, res[0], res[1][0], res[1][1], res[1][2], res[1][3])
    #print("Does the (aggregate) proof verify?", verify)

    return res


def proof_sigi0(gsp, X, Y, M1, M2, t):
    """ create GS proof that sig verifies with the signature and all but the first message secret """

    res = []

    #print("----- first equation")
    #res1 = e(R,V) * e(S,g2)* e(g1,Z)^-1 * Sum[] e(m1[i],W[i]) ]
    res1=[]
    B1 = []
    A1 = []
    C1 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            c = Bn(0)      
            row.append(c)
        C1.append(row)
    C1[0][1] = Bn(-1) # e(g1,Z)^-1
    C1[2+len(M2)][2] = Bn(1) # e(R,V)
    C1[3+len(M2)][0] = Bn(1) # e(S,g2)
    C1[0][3] = t # e(g1^t,W0)
    C1[1+len(M2)][4] = Bn(1) # e(h_i,W1)

    success1, result1 = gsp.Prove_aggreg("PPE", X, B1, A1, Y, C1, GTElem.zero(gsp.G))
    verify1 = 0
    if success1:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result1
        #verify1 = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify1:
        res1 = [C1, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res1)
    #print("Do we successfully create a first proof?", success1)
    #print("Does the first proof successfully verify?", verify1)
    
    #print("----- second equation")
    #res2 = e(R,T) * e(g1,g2)^-1 * Sum [e(U[j],m2[j])]
    res2=[]
    B2 = []
    A2 = []
    C2 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            c = Bn(0)
            row.append(c)
        C2.append(row)
    C2[0][0] = Bn(-1) # e(g1,g2)^-1
    C2[2+len(M2)][5] = Bn(1) # e(R,T)
    for i in range(len(M2)):
        C2[1+i][len(Y)-len(M2)+i] = Bn(1) # e(U, M2)
    
    success2, result2 = gsp.Prove_aggreg("PPE", X, B2, A2, Y, C2, GTElem.zero(gsp.G))
    verify2 = 0
    if success2:
        eq_type, X2, Y2, C2, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result2
        #verify2 = gsp.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify2:
        res2 = [C2, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res2)
    #print("Do we successfully create a second proof?", success2)
    #print("Does the second proof successfully verify?", verify2)
    #print("Are all the proofs successfull created?", success1*success2)
    #print("Do all the proofs verify?", verify1*verify2)

    b1 = challenge(result1)
    b2 = challenge(result2)
    C = []
    for i in range(len(C1)):
        row = []
        for j in range(len(C1[0])):
            cij = Bn(0)
            if C1[i][j] != 0:
                cij += b1 * C1[i][j]
            if C2[i][j] != 0:
                cij += b2 * C2[i][j]
            row.append(cij)
        C.append(row)
    pi = []
    for i in range(len(res1[1])):
        pi_i = []
        for j in range(len(res1[1][0])):
            pi_ij = b1*res1[1][i][j] + b2*res2[1][i][j]
            pi_i.append(pi_ij)
        pi.append(pi_i)

    #verify = gsp.Verify("PPE", X, Y, C, pi[0], pi[1], pi[2], pi[3])
    #print("Does the (aggregate) proof verify?", verify)

    return C, pi, res

def proof_sigit(gsp, X, Y, M1, M2, t):
    """ create GS proof that sig verifies with the signature and all but the first message secret """

    res = []

    #print("----- first equation")
    #res1 = e(R,V) * e(S,g2)* e(g1,Z)^-1 * Sum[] e(m1[i],W[i]) ]
    res1 = []
    B1 = []
    A1 = []
    C1 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            c = Bn(0)      
            row.append(c)
        C1.append(row)
    C1[0][1] = Bn(-1) # e(g1,Z)^-1
    C1[4+len(M2)][2] = Bn(1) # e(R,V)
    C1[5+len(M2)][0] = Bn(1) # e(S,g2)
    C1[0][3] = t # e(g1^t,W0)
    C1[1+len(M2)][4] = Bn(1) # e(h_i,W1)

    success1, result1 = gsp.Prove_aggreg("PPE", X, B1, A1, Y, C1, GTElem.zero(gsp.G))
    verify1 = 0
    if success1:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result1
        #verify1 = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify1:
        res1 = [C1, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res1)
    #print("Do we successfully create a first proof?", success1)
    #print("Does the first proof successfully verify?", verify1)
    
    #print("----- second equation")
    #res2 = e(R,T) * e(g1,g2)^-1 * Sum [e(U[j],m2[j])]
    res2 = []
    B2 = []
    A2 = []
    C2 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            c = Bn(0)
            row.append(c)
        C2.append(row)
    C2[0][0] = Bn(-1) # e(g1,g2)^-1
    C2[4+len(M2)][6] = Bn(1) # e(R,T)
    for i in range(len(M2)):
        C2[1+i][len(Y)-len(M2)+i] = Bn(1) # e(U, M2)
    
    success2, result2 = gsp.Prove_aggreg("PPE", X, B2, A2, Y, C2, GTElem.zero(gsp.G))
    verify2 = 0
    if success2:
        eq_type, X2, Y2, C2, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result2
        #verify2 = gsp.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify2:
        res2 = [C2, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res2)
    #print("Do we successfully create a second proof?", success2)
    #print("Does the second proof successfully verify?", verify2)

    #print("Are all the proofs successfull created?", success1*success2)
    #print("Do all the proofs verify?", verify1*verify2)

    b1 = challenge(result1)
    b2 = challenge(result2)
    C = []
    for i in range(len(C1)):
        row = []
        for j in range(len(C1[0])):
            cij = Bn(0)
            if C1[i][j] != 0:
                cij += b1*C1[i][j]
            if C2[i][j] != 0:
                cij += b2*C2[i][j]
            row.append(cij)
        C.append(row)
    pi = []
    for i in range(len(res1[1])):
        pi_i = []
        for j in range(len(res1[1][0])):
            pi_ij = b1*res1[1][i][j] + b2*res2[1][i][j]
            pi_i.append(pi_ij)
        pi.append(pi_i)

    #verify = gsp.Verify("PPE", X, Y, C, pi[0], pi[1], pi[2], pi[3])
    #print("Does the (aggregate) proof verify?", verify)

    return C, pi, res

def enc_proof(gsp, X, Y, a):
    res = []

    #print("--- first equation")
    # e = h_i + h*r
    res1 = []
    A1 = []
    B1 = []
    C1 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C1.append(row)
    C1[12][0] = Bn(-1) # - e
    C1[7][0] = Bn(1) # + h_i
    C1[4][2] = Bn(1) # + h*r

    success1, result1 = gsp.Prove_aggreg("ME1", X, B1, A1, Y, C1, G1Elem.inf(gsp.G))
    verify1 = 0
    if success1:
        eq_type, X1, Y1, C1, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result1
        #verify1 = gsp.Verify(eq_type, X1, Y1, C1, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify1:
        res1 = [C1, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res1)
    #print("Do we successfully create a first proof?", success1)
    #print("Does the first proof successfully verify?", verify1)

    #print("--- second equation")
    #u1 = g1_enc*r
    res2 = []
    B2 = []
    A2 = []
    C2 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C2.append(row)
    C2[0][2] = Bn(-1) # - g1enc*r
    C2[10][0] = Bn(1) # + u1
    success2, result2 = gsp.Prove_aggreg("MC1", X, B2, A2, Y, C2, G1Elem.inf(gsp.G))
    verify2 = 0
    if success2:
        eq_type, X2, Y2, C2, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result2
        #verify2 = gsp.Verify(eq_type, X2, Y2, C2, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify2:
        res2 = [C2, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res2)
    #print("Do we successfully create a second proof?", success2)
    #print("Does the second proof successfully verify?", verify2)

    #print("--- third equation")
    #u2 = g2_enc*r
    res3 = []
    B3 = []
    A3 = []
    C3 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C3.append(row)
    C3[1][2] = Bn(-1) # - g2_enc * r
    C3[11][0] = Bn(1) # + u2
    success3, result3 = gsp.Prove_aggreg("ME1", X, B3, A3, Y, C3, G1Elem.inf(gsp.G))
    verify3 = 0
    if success3:
        eq_type, X3, Y3, C3, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result3
        #verify3 = gsp.Verify(eq_type, X3, Y3, C3, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify3:
        res3 = [C3, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res3)
    #print("Do we successfully create a third proof?", success3)
    #print("Does the third proof successfully verify?", verify3)

    """ We perform this check with a hash
    #print("--- fourth equation")
    # da = d*a
    res4 = []
    B4 = []
    A4 = []
    C4 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C4.append(row)
    C4[9][0] = Bn(-1) # - da
    C4[3][0] = a # + d*a
    success4, result4 = gsp.Prove_aggreg("ME1", X, B4, A4, Y, C4, G1Elem.inf(gsp.G))
    verify4 = 0
    if success4:
        eq_type, X4, Y4, C4, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result4
        #verify4 = gsp.Verify(eq_type, X4, Y4, C4, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify4:
        res4 = [C4, [pi2_v1, pi2_w1, pi1_v2, pi1_w2] ]
        res.append(res4)
    #print("Do we successfully create an fourth proof?", success4)
    #print("Does the fourth proof successfully verify?", verify4)
    """

    #print("--- fifth equation")
    #v = c*r + d*(a*r) = c*r + da*r
    res5 = []
    B5 = []
    A5 = []
    C5 = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C5.append(row)
    C5[2][2] = Bn(1) # + c*r
    C5[9][2] = Bn(1) # + da*r
    C5[13][0] = Bn(-1) # - v
    success5, result5 = gsp.Prove_aggreg("MC1", X, B5, A5, Y, C5, G1Elem.inf(gsp.G))
    verify5 = 0
    if success5:
        eq_type, X5, Y5, C5, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result5
        #verify5 = gsp.Verify(eq_type, X5, Y5, C5, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify5:
        res5 = [C5, [pi2_v1, pi2_w1, pi1_v2, pi1_w2]]
        res.append(res5)
    #print("Do we successfully create a fifth proof?", success5)
    #print("Does the fifth proof successfully verify?", verify5)

    #print("Do we successfully create all the proofs?", success1*success2*success3*success4*success5)
    #print("Do all the proofs successfully verify?", verify1*verify2*verify3*verify4*verify5)

    b1 = challenge(result1)
    b2 = challenge(result2)
    b3 = challenge(result3)
    #b4 = challenge(result4)
    b5 = challenge(result5)
    C = []
    for i in range(len(C1)):
        row = []
        for j in range(len(C1[0])):
            cij = Bn(0)
            if C1[i][j] != 0:
                cij += b1*C1[i][j]
            if C2[i][j] != 0:
                cij += b2*C2[i][j]
            if C3[i][j] != 0:
                cij += b3*C3[i][j]
            #if C4[i][j] != 0:
            #    cij += b4*C4[i][j]
            if C5[i][j] != 0:
                cij += b5*C5[i][j]
            row.append(cij)
        C.append(row)

    pi = []
    for i in range(len(res1[1])):
        pi_i = []
        for j in range(len(res1[1][0])):
            pi_ij = b1*res1[1][i][j] + b2*res2[1][i][j] + b3*res3[1][i][j] + b5*res5[1][i][j] # + b4*res4[1][i][j]
            pi_i.append(pi_ij)
        pi.append(pi_i)

    print("\n--- enc")
    for i in range(len(res)):
        print("proof #"+str(i))
        for j in range(4):
            for k in range(2):
               if type(res[i][1][j][k]) == G1Elem:
                   print(i,j,k, type(res[i][1][j][k]), res[i][1][j][k].eq(G1Elem.inf(gsp.G)))
               else:
                   print(i,j,k,type(res[i][1][j][k]), res[i][1][j][k].eq(G2Elem.inf(gsp.G)))

    #verify = gsp.Verify("ME1", X, Y, C, pi[0], pi[1], pi[2], pi[3])
    #print("Does the (aggregate) proof verify?", verify)

    return C, pi, res

def proof_pkuv(gsp, X, Y):
    res = []
    B = []
    A = []
    C = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C.append(row)
    C[6][0] = Bn(1) # pk_uv
    C[5][1] = Bn(-1) # - h_v^sku
    success, result = gsp.Prove_aggreg("MC1", X, B, A, Y, C, G1Elem.inf(gsp.G))
    verify = 0
    if success:
        eq_type, X, Y, C, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result
        #verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify:
        res = [C, [pi2_v1, pi2_w1, pi1_v2, pi1_w2] ]
    #print("Do we successfully create a proof?", success)
    #print("Does the proof successfully verify?", verify)

    print("\n----pk uv")
    for i in range(len(res[1])):
        for j in range(2):
           if type(res[1][i][j]) == G1Elem:
               print(i,j,type(res[1][i][j]), res[1][i][j].eq(G1Elem.inf(gsp.G)))
           else:
               print(i,j,type(res[1][i][j]), res[1][i][j].eq(G2Elem.inf(gsp.G)))

    b = challenge(result)
    for i in range(len(res[0])):
        for j in range(len(res[0][0])):
            if res[0][i][j] != 0:
                res[0][i][j]= b*res[0][i][j]
    for i in range(len(res[1])):
        for j in range(len(res[1][i])):
            res[1][i][j] = b*res[1][i][j]

    #verify = gsp.Verify("ME1", X, Y, C, res[1][0], res[1][1], res[1][2], res[1][3])
    #print("Does the (aggregate) proof verify?", verify)

    return res

def proof_pkui(gsp, X, Y):
    res = []
    B = []
    A = []
    C = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):
            row.append(Bn(0))
        C.append(row)
    C[8][0] = Bn(1) # pk_ui
    C[7][1] = Bn(-1) # - h_i^sku
    success, result = gsp.Prove_aggreg("ME1", X, B, A, Y, C, G1Elem.inf(gsp.G))
    verify = 0
    if success:
        eq_type, X, Y, C, T_eq, pi2_v1, pi2_w1, pi1_v2, pi1_w2 = result
        #verify = gsp.Verify(eq_type, X, Y, C, pi2_v1, pi2_w1, pi1_v2, pi1_w2)
        #if verify:
        res = [C, [pi2_v1, pi2_w1, pi1_v2, pi1_w2] ]
    #print("Do we successfully create a proof?", success)
    #print("Does the proof successfully verify?", verify)

    print("\n--- pkui")
    for i in range(len(res[1])):
        for j in range(2):
           if type(res[1][i][j]) == G1Elem:
               print(i,j,type(res[1][i][j]),res[1][i][j].eq(G1Elem.inf(gsp.G)))
           else:
               print(i,j,type(res[1][i][j]),res[1][i][j].eq(G2Elem.inf(gsp.G)))

    b = challenge(result)
    for i in range(len(res[0])):
        for j in range(len(res[0][0])):
            res[0][i][j]= b*res[0][i][j]
    for i in range(len(res[1])):
        for j in range(len(res[1][i])):
            res[1][i][j] = b*res[1][i][j]

    #verify = gsp.Verify("ME1", X, Y, res[0], res[1][0], res[1][1], res[1][2], res[1][3])
    #print("Does the (aggregate) proof verify?", verify)

    return res

def challenge(elements):
    """Packages a challenge in a bijective way"""
    elem = [len(elements)] + elements
    elem_str = map(str, elem)
    elem_len = map(lambda x: "%s||%s" % (len(x) , x), elem_str)
    state = "|".join(elem_len)
    H = sha256()
    H.update(state.encode("utf8"))
    return Bn.from_binary(H.digest())

def prepare_proofs(auth, vkI, pki, pkv, m, sig_t, t, pk_ui, pk_uv, sku, cipher, ek, r):
    #print("Prepare aggregated proofs")
    #print("Prepare proof: Commit")
    #import time
    #t_start = time.time()
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
    cm_sig0_i.append(auth.GS.Commit({"group":1, "type":"com", "value":sig0_i[0]}))
    cm_sig0_i.append(auth.GS.Commit({"group":1, "type":"enc", "value":sig0_i[1]}))
    cm_sig0_i.append(auth.GS.Commit({"group":2, "type":"enc", "value":sig0_i[2]}))

    h_v, vk_v, sig0_v = pkv
    cm_h_v = auth.GS.Commit({"group":1, "type":"pub", "value":h_v})

    cm_m = []
    for i in range(len(m)):
        cm_m.append(auth.GS.Commit({"group":1, "type":"pub", "value":m[i]}))
    for i in range(len(cm_vk_i) - 1 -2 - len(cm_m)):
        cm_m.append(auth.GS.Commit({"group":1, "type":"pub", "value":G1Elem.inf(auth.GS.G)}))

    sig_i, c_it = sig_t
    cm_sig_i = []
    for i in range(len(sig_i)):
        cm_sig_i.append(auth.GS.Commit({"group":1, "type":"enc", "value":sig_i[i]}))
    _, sigt_i = c_it
    cm_sigt_i = []
    cm_sigt_i.append(auth.GS.Commit({"group":1, "type":"com", "value":sigt_i[0]}))
    cm_sigt_i.append(auth.GS.Commit({"group":1, "type":"enc", "value":sigt_i[1]}))
    cm_sigt_i.append(auth.GS.Commit({"group":2, "type":"enc", "value":sigt_i[2]}))

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

    #t_commit = time.time()
    #print("--- Commitment time:", t_commit - t_start)

    x_ppe = [cm_g1]
    x_ppe.extend(cm_vkI[0])
    x_ppe.append(cm_h_i)
    x_ppe.extend([cm_sig0_i[0], cm_sig0_i[1]])
    x_ppe.extend([cm_sigt_i[0], cm_sigt_i[1]])
    x_ppe.extend(cm_sig_i)
    x_ppe.append(cm_pk_ui)
    x_ppe.extend(cm_m)
    #print("\nx_ppe", len(x_ppe), x_ppe)

    y_ppe = [cm_g2]
    y_ppe.append(cm_vkI[3])
    y_ppe.append(cm_vkI[2])
    y_ppe.extend(cm_vkI[1])
    y_ppe.append(cm_sig0_i[2])
    y_ppe.append(cm_sigt_i[2])
    y_ppe.extend(cm_vk_i)
    #print("\ny_ppe", len(y_ppe), y_ppe)

    x_me1 = []
    x_me1.extend(cm_params_enc)
    x_me1.extend(cm_ek)
    x_me1.append(cm_h_v)
    x_me1.append(cm_pk_uv)
    x_me1.append(cm_h_i)
    x_me1.append(cm_pk_ui)
    x_me1.append(cm_da)
    x_me1.extend(cm_cipher)
    #print("\nx_me1", len(x_me1), x_me1)

    y_me1 = [cm_1, cm_sku, cm_r[0]]
    #print("\ny_me1", len(y_me1), y_me1)

    #print("Prepare proof: Prove")
    #t_prove = time.time()
    #print("--- sigi proof")
    cm_msg = [cm_pk_ui]
    cm_msg.extend(cm_m)
    for i in range(len(cm_vk_i) - 2 - len(cm_msg)):
        cm_msg.append(auth.GS.Commit({"group":1, "type":"pub", "value":G1Elem.inf(auth.GS.G)}))
    C_sigi, pi_sigi = proof_sigi(auth.GS, x_ppe, y_ppe, cm_msg)

    cm_msg1 = [cm_g1, cm_h_i]
    cm_msg2 = cm_vk_i

    #print("--- sigi0 proof")
    C_sigi0, pi_sigi0, _ = proof_sigi0(auth.GS, x_ppe, y_ppe, cm_msg1, cm_msg2, Bn(0))

    #print("--- sigit proof")
    C_sigit, pi_sigit, _ = proof_sigit(auth.GS, x_ppe, y_ppe, cm_msg1, cm_msg2, t)

    #print("--- Aggregate PPE proofs")
    c_ppe = []
    for i in range(len(C_sigi)):
        row = []
        for j in range(len(C_sigi[i])):
            cij = C_sigi[i][j] + C_sigit[i][j] + C_sigi0[i][j]
            row.append(cij)
        c_ppe.append(row)

    pi_ppe = []
    for i in range(len(pi_sigi)):
        pi_i = []
        for j in range(len(pi_sigi[i])):
            pi_ij = pi_sigi[i][j] + pi_sigit[i][j] + pi_sigi0[i][j]
            pi_i.append(pi_ij)
        pi_ppe.append(pi_i)
    
    #print("--- Randomize PPE proof")
    pi_ppe[0], pi_ppe[1], pi_ppe[2], pi_ppe[3] = auth.GS.Randomize("PPE", pi_ppe[0], pi_ppe[1], pi_ppe[2], pi_ppe[3])
    res_ppe = [c_ppe, pi_ppe]
    #t_ppe = time.time()

    print("--- exponent proofs pk_ui")
    C_pkui, pi_pkui = proof_pkui(auth.GS, x_me1, y_me1)

    print("--- exponent proofs pk_uv")
    C_pkuv, pi_pkuv = proof_pkuv(auth.GS, x_me1, y_me1)

    #print("--- enc proof")
    C_enc, pi_enc, _ = enc_proof(auth.GS, x_me1, y_me1, r[1])

    #print("--- aggregate ME1 proofs")
    c_me1 = []
    for i in range(len(C_enc)):
        row = []
        for j in range(len(C_enc[i])):
            cij = C_enc[i][j] + C_pkui[i][j] + C_pkuv[i][j]
            row.append(cij)
        c_me1.append(row)

    pi_me1 = []
    for i in range(len(pi_enc)):
        pi_i = []
        for j in range(len(pi_enc[i])):
            pi_ij = pi_enc[i][j] + pi_pkui[i][j] + pi_pkuv[i][j]
            pi_i.append(pi_ij)
        pi_me1.append(pi_i)

    #print("------ Randomize ME1 proof")
    pi_me1[0], pi_me1[1], pi_me1[2], pi_me1[3] = auth.GS.Randomize("PPE", pi_me1[0], pi_me1[1], pi_me1[2], pi_me1[3])
    res_me1 = [c_me1, pi_me1]

    #t_end = time.time()
    #print("--- Prove & aggregation time:", t_end - t_commit, "(PPE proof: "+ str(t_ppe-t_commit) +"+ ME1 proof"+ str(t_end-t_ppe) +")")

    verify_ppe = auth.GS.Verify(
        "PPE",
        x_ppe,
        y_ppe,
        c_ppe,
        pi_ppe[0],
        pi_ppe[1],
        pi_ppe[2],
        pi_ppe[3]
        )
    print("------ Aggregate PPE verify?", verify_ppe)

    verify_me1 = auth.GS.Verify(
        "ME1",
        x_me1,
        y_me1,
        c_me1,
        pi_me1[0],
        pi_me1[1],
        pi_me1[2],
        pi_me1[3]
        )
    print("------ Aggregate ME1 verify?", verify_me1)

    verify = verify_ppe*verify_me1
    res = [
        [res_ppe[1], ["PPE", x_ppe, y_ppe, res_ppe[0], GTElem.zero(auth.GS.G)]],
        [res_me1[1], ["ME1", x_me1, y_me1, res_me1[0], G1Elem.inf(auth.GS.G)]]
        ]
    #print("--- Do all the proofs verify?", verify_ppe*verify_me1)

    return verify, res
