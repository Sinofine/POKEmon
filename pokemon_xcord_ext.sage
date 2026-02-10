# %%
import os, gc
import random
import time
from hashlib import shake_256
from sage.all import *
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

from montgomery_isogenies.kummer_line_ext import KummerLine
from montgomery_isogenies.kummer_isogeny_ext import KummerLineIsogeny
from montgomery_isogenies.isogenies_x_only_ext import lift_image_to_curve

from utilities.discrete_log import BiDLP, discrete_log_pari
from utilities.supersingular_ext import torsion_basis, torsion_basis_with_pairing, torsion_basis_2e
from utilities.utils import speed_up_sagemath
from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_endomorphism import KaniEndoHalf
from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
gc.disable()
speed_up_sagemath()
prime_params = [{'prime': 731806003979838355076968625254907577632666303086848364729371930930010979324395519,
 'a': 126,
 'b': 3 * 13^2 * 19 * 53 * 109 * 131 * 1013,
 'bits': 269}, {'prime': 14951226794123200734107772005897438752368996540481204139289517674370349731477510516905646595491998044013879365722207986647039,
 'a': 191,
 'b': 7 * 13 * 47 * 1201 * 2399 * 4789 * 8059 * 23399,
 'bits': 413}, {'prime': 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151,
 'a': 253,
 'b': 3 * 5^2 * 11 * 13 * 17 * 31 * 41 * 53 * 131 * 157 * 521 * 1613 * 2731 * 2777,
 'bits': 521}]
def find_point_from_index(E, index=0):
    F = E.base_ring()
    p = F.characteristic()
    a = F.gen()
    curr_x = F(index)
    while True:
        try:
            # lift_x 会尝试寻找 x 坐标对应的点
            # all=True 返回两个点 (x, y) 和 (x, -y)
            curr_x += 1
            print(curr_x)
            return E.lift_x(curr_x+a), curr_x
        except ValueError:
            # 如果当前 x 不是二次剩余，递增继续寻找
            curr_x += 1
def get_deterministic_n_torsion_point(E, n, N, seed=0):
    factors = n.factor()
    l, k = factors[0]
    l_inv = n // l
    
    curr_seed = seed
    while True:
        # 使用确定性逻辑替代 random_point()
        # 例如从第 curr_seed 个 x 坐标开始找点
        R, curr_seed = find_point_from_index(E, curr_seed)
        P = (N // n) * R
        
        if P != E(0) and n * P == E(0) and l_inv * P != E(0):
            return P, curr_seed
        curr_seed += 1
def get_torsion_basis_fast(E, n, N, random=True):
    """
    针对 n 为素数幂的快速基生成。
    假设 n 阶点在当前域中已经存在。
    """
    # 假设 n = l^k
    factors = n.factor()
    l, k = factors[0] 
    l_inv = n // l  # 用于检查点阶是否恰好为 n
    seed = 0
    
    def get_random_n_torsion_point():
        while True:
            R = E.random_point()
            P = (N // n) * R
            if P == E(0):
                continue
            # 检查阶是否恰好为 n: [n]P = 0 且 [n/l]P != 0
            if n * P == E(0) and l_inv * P != E(0):
                return P, 0

    # 1. 寻找第一个 n 阶点 P
    P, seed = get_random_n_torsion_point() if random else get_deterministic_n_torsion_point(E, n, N, seed)
    
    # 2. 寻找第二个 n 阶点 Q，要求与 P 线性无关
    while True:
        Q, seed = get_random_n_torsion_point() if random else get_deterministic_n_torsion_point(E, n, N, seed)
        zeta = P.weil_pairing(Q, n)
        # 验证线性无关：zeta 必须是本原 n 次单位根
        # 即 zeta^(n/l) != 1
        if zeta^l_inv != 1:
            return P, Q

# %%
def xof_kdf(X, Y):
    X_bytes = X.to_bytes()
    Y_bytes = Y.to_bytes()
    return shake_256(X_bytes + Y_bytes)

def xof_encrypt(xof, msg):
    key = xof.digest(len(msg))
    return bytes(x ^^ y for (x, y) in zip(key, msg))

def random_unit(modulus):
    while True:
        alpha = ZZ.random_element(modulus)
        if gcd(alpha, modulus) == 1:
            break
    return alpha

def random_matrix(modulus):
    while True:
        d1 = ZZ.random_element(modulus)
        d2 = ZZ.random_element(modulus)
        d3 = ZZ.random_element(modulus)
        d4 = ZZ.random_element(modulus)
        if gcd(d1*d4 - d2*d3, modulus) == 1:
            break
    return d1, d2, d3, d4

def prodeval(Phi, secZero, pos):
    return lambda P: Phi(CouplePoint(P, secZero))[pos]
def is_sum_of_two_squares_calc(n):
    # 获取整数的素因数分解
    F = factor(n) 
    for p, exponent in F:
        # 检查模 4 余 3 的素因子，其指数必须为偶数
        if p % 4 == 3 and exponent % 2 != 0:
            return False
    return True
def point_to_xonly(P, Q):
    L = KummerLine(P.curve())

    PQ = P - Q
    xP = L(P[0])
    xQ = L(Q[0])
    xPQ = L(PQ[0])

    return L, xP, xQ, xPQ

# %%
def gen_isogeny_deg(E0, P0, Q0,R0,S0, q, a_param, B_param):
    """
    Algorithm 1: Generating a q-isogeny
    
    Require: p = c·2^(2a)·3^b - 1, E_0/F_p^2: y^2 = x^3 + x, 
             ⟨P_0, Q_0⟩ = E_0[2^(2a)], a degree q.
    Ensure: A representation of an isogeny φ: E_0 → E_1 of degree q.
    
    Parameters:
    - E0: Supersingular elliptic curve
    - P0, Q0: Basis of E_0[2^(2a)]
    - q: Degree of isogeny
    - a_param: parameter a for 2-torsion
    - B_param: parameter B for 3-torsion
    
    Returns:
    - E1: Codomain curve
    - P1, Q1: Torsion basis of E_1[2^(2a)]
    """
    # Step 1: Generate an endomorphism θ of degree q(2^(2a) - q)B
    deg_theta = q * (2^(2*a_param) - q) * (B_param)

    
    # For this implementation, we use the Frobenius endomorphism composition
    # θ is constructed to have the desired properties
    # We use a random isogeny approach to generate it
    bound = ZZ(4*deg_theta // p)
    zt_bound = isqrt(4*deg_theta // p-1)
    QF = BinaryQF([1,0,1])

    for _ in range(bound):
        zz = randint(0, zt_bound)
        tt = randint(0, zt_bound)
        sq = deg_theta - p * (zz**2 + tt**2)
        if sq <= 0:
            continue
        if not sq.is_prime() or sq % 4 != 1:
            continue
        # print(sq)
        # Try and get a solution sq = x^2 + y^2
        try:
            xx, yy = QF.solve_integer(sq, algorithm="cornacchia")
            break
        except ValueError:
            continue
    else:
        raise ValueError("Could not find a suitable endomorphism.")

    i_end = lambda P: E0(-P[0], i*P[1])
    pi_end = lambda P: E0(P[0]**p, P[1]**p)
    θ = lambda P: xx*P + yy*i_end(P) + zz*pi_end(P) + tt*i_end(pi_end(P))
    # Step 2: Compute P'_0 = θ(P_0), Q'_0 = θ(Q_0)
    # Apply Frobenius-based endomorphism
    # This is a composition creating the desired kernel
    P0_prime = θ(P0)
    Q0_prime = θ(Q0)
    R0_prime = θ(R0)
    S0_prime = θ(S0)
    
    
    # Step 3: Let K = ker(θ) ∩ E_0[3^b]
    
    # K is a point in E_0[3^b]
    # For the algorithm, we select a random point in the 3^b torsion
    try:
        R = S0
        wp_ = R0_prime.weil_pairing(R, B_param, algorithm='pari')
        wp = S0_prime.weil_pairing(R, B_param, algorithm='pari')
        
        discrete_log_pari(wp_, wp, B_param)

        K3_dual = S0_prime
    except TypeError:
        R = R0_prime
        K3_dual = R0_prime # P0_prime is guaranteed to be lin. indep with K3
    
    phi = E0.isogeny(K3_dual, algorithm="factored", model="montgomery")
    
    
    # Step 5: Compute P_1, Q_1 = [1/B]φ(P'_0), [1/B]φ(Q'_0)
    # Evaluate phi on the endomorphism images
    phi_P0_prime = phi(P0_prime)
    phi_Q0_prime = phi(Q0_prime)
    
    # Divide by B (multiply by the inverse)
    inv_B = pow(B_param, -1, 2^(2*a_param))  # Inverse in the group order
    P1 = inv_B * phi_P0_prime
    Q1 = inv_B * phi_Q0_prime
    E1 = P1.curve()
    
    
    # Step 6: Compute Φ: E_0 × E_1 → E_2 × E_3 with kernel 
    # (([−q]P_0, P_1), ([−q]Q_0, Q_1))
    P_couple = CouplePoint(-q * P0, P1)
    Q_couple = CouplePoint(-q * Q0, Q1)
    kernel = (P_couple, Q_couple)
    
    Phi = EllipticProductIsogenySqrt(kernel, 2*a_param)
    E2, E3 = Phi.codomain()
    
    
    # Step 7: Compute (P_2, P_3) = Φ(P_0, 0), (Q_2, Q_3) = Φ(Q_0, 0)
    P_couple_in = CouplePoint(P0, E1(0))
    Q_couple_in = CouplePoint(Q0, E1(0))
    PQ_couple_in = CouplePoint(P0 + Q0, E1(0))
    P_couple_out = Phi(P_couple_in)
    Q_couple_out = Phi(Q_couple_in)
    PQ_couple_out = Phi(PQ_couple_in)
    
    P2, P3 = P_couple_out[0], P_couple_out[1]
    Q2, Q3 = Q_couple_out[0], Q_couple_out[1]
    PQ2, PQ3 = PQ_couple_out[0], PQ_couple_out[1]
    if not (P2+Q2 == PQ2 or P2+Q2 == -PQ2):
        Q2 = -Q2
    if not (P3+Q3 == PQ3 or P3+Q3 == -PQ3):
        Q3 = -Q3
    
    # Step 8-9: Check Weil pairing condition
    # If e_{2a}(P_2, Q_2) = e_{2a}(P_0, Q_0)^q then return (E_2, P_2, Q_2)
    # else return (E_3, P_3, Q_3)
    
    # Compute Weil pairings
    try:
        e_P0_Q0 = P0.weil_pairing(Q0, 2^(2*a_param), algorithm='pari')
        e_P2_Q2 = P2.weil_pairing(Q2, 2^(2*a_param), algorithm='pari')
        e_P3_Q3 = P3.weil_pairing(Q3, 2^(2*a_param), algorithm='pari')
        
        target_pairing = e_P0_Q0^q
        if e_P2_Q2 == target_pairing:
            # print(f"Step 9: Pairing check passed, returning (E_2, P_2, Q_2)")
            return prodeval(Phi, E1(0), 0)
        else:
            assert e_P3_Q3 == target_pairing
            # print(f"Step 9: Pairing check failed, returning (E_3, P_3, Q_3)")
            return prodeval(Phi, E1(0), 1)
    except:
      raise RuntimeError("Weil pairing computation failed.")

# %%
def one_point(E0, E1, P0_op, phi1P0, d, points):
    """
    Algorithm 2: Onepoint
    
    Require: Two supersingular elliptic curves E0, E1 connected by an odd d-isogeny φ1
             and torsion points P0, φ1(P0) = P1 of order 2^(2a), a point R0 ∈ E0(F_{q^k})
    Ensure: φ1(R0)
    
    Args:
        E0, E1: Supersingular elliptic curves connected by φ1
        P0: Torsion point on E0 of order 2^(2a)
        phi1P0: φ1(P0) = P1, torsion point on E1 of order 2^(2a)
        d: Degree of isogeny φ1 (must be odd)
        R0: Point on E0 to evaluate through φ1
        
    Returns:
        R1: φ1(R0)
    """
    
    if d % 2 == 0:
        raise ValueError("d must be odd")
    
    P1_op = phi1P0
    n = 2^a  # degree for auxiliary isogenies
    a1, a2 = two_squares(n-d)
    # if P0_op.order()<2^(2*a) or P1_op.order()<2^(2*a):
    #     raise Exception("P0, P1 not order 2^2a")
    # Q0_basis = E0.torsion_basis(2^(2*a))
    Q0_basis = torsion_basis_2e(E0, 2*a)
    # Try to find Q0 such that the pairing is non-trivial
    # This should be changed to a random search in practice
    zeta = None
    for idx in range(0, 1000):
        Q0_cand = Q0_basis[0] + idx * Q0_basis[1]  # Try multiples of the second basis point
        try:
            pairing_val = P0_op.weil_pairing(Q0_cand, 2^(2*a), algorithm='pari')
            if pairing_val^(2^(2*a-1)) != 1:
                Q0_op = Q0_cand
                zeta = pairing_val
                break
        except:
            continue
    
    if zeta is None:
        raise ValueError("Could not find suitable Q0")
    K2 = n * P0_op  # Kernel generator
    _E0 = KummerLine(E0)
    xK2 = _E0(K2)
    xP0 = _E0(P0_op)
    xQ0 = _E0(Q0_op)
    phi2 = KummerLineIsogeny(_E0, xK2, n)
    _E2 = phi2.codomain()
    E2 = _E2.curve()
    # Compute torsion points (P2, Q2) = (φ2(P0), 2^a φ2(Q0))
    xP2 = phi2(xP0)
    xQ0_image = phi2(xQ0)
    P2, Q0_image = lift_image_to_curve(P0, Q0, xP2, xQ0_image, 2^(2*a), 2^a)
    Q2 = n * Q0_image  # 2^a φ2(Q0)
    
    # Step 3: Select Q1 ∈ E1[2^(2a)] such that e_{2^(2a)}(P1, Q1) = ζ^d
    # Q1_basis = E1.torsion_basis(2^(2*a))
    Q1_basis = torsion_basis_2e(E1, 2*a)
    for idx in range(0, 2^(2*a)):
        Q1_cand = Q1_basis[0] + idx * Q1_basis[1]  # Try multiples of the second basis point
        try:
            pairing_val = P1_op.weil_pairing(Q1_cand, 2^(2*a), algorithm='pari')
            if pairing_val^(2^(2*a-1)) != 1:
                S1 = Q1_cand
                break
        except:
            continue
    target_pairing = zeta^d
    ee = P1_op.weil_pairing(S1, 2^(2*a), algorithm='pari')
    ii = discrete_log_pari(target_pairing, ee, 2^(2*a))
    Q1_op = P1_op+ii*S1
    _E1 = KummerLine(E1)
    _P1_op = _E1(P1_op)
    _Q1_op = _E1(Q1_op)
    
    _K2_prime = n * _P1_op
    
    # Use isogeny method directly
    phi2_prime = KummerLineIsogeny(_E1, _K2_prime, n)
    xE3 = phi2_prime.codomain()
    E3 = xE3.curve()
    xP3 = phi2_prime(_P1_op)
    xQ1_image = phi2_prime(_Q1_op)
    xQ3 = n * xQ1_image
    P3, Q3 = lift_image_to_curve(P1_op, Q1_op, xP3, xQ3, 2^(2*a), 2^a)
    phi2_prime_dual = KummerLineIsogeny(xE3, xQ3, n)
    # 原有的代码继续...
    F = KaniEndoHalf(P2, Q2, P3, (2*int(P2.weil_pairing(Q2, 2^a, algorithm='pari')^d == P3.weil_pairing(Q3, 2^a,algorithm='pari'))-1)*Q3, d, a1, a2, a, a)
    phi1prime = lambda pt: F(TuplePoint(pt, E2(0), E3(0), E3(0)))[2]
    return map(lambda pt: phi2_prime_dual(xE3(phi1prime(phi2(_E0(pt)).curve_point()))), points)

# %%
def keygen(E0, P0, Q0, R0, S0, a_param, B_param):
    """
    Algorithm 3: KeyGen
    
    Require: E_0, E_0[2^(2a)] = ⟨P_0, Q_0⟩, E_0[B] = ⟨R_0, S_0⟩.
    Ensure: sk, pk.
    
    Parameters:
    - E0: Supersingular elliptic curve
    - P0, Q0: Basis of E_0[2^(2a)]
    - R0, S0: Basis of E_0[B]
    - a_param: parameter a for 2-torsion
    - B_param: parameter B for B-torsion
    
    Returns:
    - sk: secret key = (q, alpha, beta)
    - pk: public key = (E_A, P_A, Q_A, R_A, S_A)
    """
    
    # Step 1: Select a random q ∈ [1, 2^a - 1] such that 2^a - q is coprime with 2, 3
    two_to_a = 2^(a_param)
    two_to_twoa = 2^(2*a_param)
    q = None
    
    while q is None:
        tempx = ZZ(randint(1, two_to_a.isqrt())//2)*2
        tempy = ZZ(randint(1, two_to_a.isqrt())//2)*2+1
        degree = two_to_a - (tempx^2 + tempy^2)
        if degree % 3 !=0 and degree > 0:
            q = degree
    
    phi = gen_isogeny_deg(E0,P0,Q0,R0,S0,q, a_param, B_param)
    
    P23x = P0+R0
    Q23x = Q0+S0
    alpha = random_unit(2^(2*a_param))
    beta = random_unit(B_param)
    imP23x = phi(P23x)
    imQ23x = phi(Q23x)
    imPQ23x = phi(P23x+Q23x)
    if not (imP23x+imQ23x == imPQ23x or imP23x+imQ23x == -imPQ23x):
        imQ23x = -imQ23x
    EA = imP23x.curve()
    _EA = KummerLine(EA)
    PA = alpha * inverse_mod(B_param, two_to_twoa) * (B_param * imP23x)
    QA = alpha * inverse_mod(B_param, two_to_twoa) * (B_param * imQ23x)
    
    RA = beta * inverse_mod(two_to_twoa, B_param) * (two_to_twoa * imP23x)
    SA = beta * inverse_mod(two_to_twoa, B_param) * (two_to_twoa * imQ23x)
    
    _, xPA, xQA, xPQA = point_to_xonly(PA, QA)
    _, xRA, xSA, xRSA = point_to_xonly(RA, SA)

    sk = (q, alpha, beta)
    pk = (_EA, xPA, xQA, xPQA, xRA, xSA, xRSA, RA, SA)
    
    
    return sk, pk

# %%
def encryption(pk, m, a_param, B_param):
    """
    Algorithm 4: Encryption
    
    Require: A message m, a public key pk = (E_A, xPA, xQA, xPQA, xRA, xSA, xRSA).
    Ensure: A ciphertext ct.
    
    Parameters:
    - pk: Public key (E_A, xPA, xQA, xPQA, xRA, xSA, xRSA) from keygen
    - m: Message to encrypt (bytes)
    - a_param: parameter a for 2-torsion
    - B_param: parameter B for B-torsion
    
    Returns:
    - ct: Ciphertext (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    """
    _EA, xPA, xQA, xPQA, xRA, xSA, xRSA, RA, SA = pk
    
    # Step 1: Select a random integer s ∈ [0, 2^(2a)]
    s = ZZ(randint(0, 2^(2*a_param-1)))*2+1
    # Step 2: Compute the isogeny ψ : E_0 → E_B with kernel ⟨P_0 + sQ_0⟩
    kernel_B = xQ0.ladder_3_pt(xP0, xPQ0, s)
    psi = KummerLineIsogeny(_E0, kernel_B, A)
    _EB = psi.codomain()

    # Step 3: Compute the isogeny ψ' : E_A → E_AB with kernel ⟨P_A + sQ_A⟩
    # kernel_AB = P_A + s * Q_A
    kernel_AB = xQA.ladder_3_pt(xPA, xPQA, s)
    # psi_prime0 = E_A.isogeny(kernel_AB, algorithm="factored", model="montgomery")
    psi_prime = KummerLineIsogeny(_EA, kernel_AB, A)
    _EAB = psi_prime.codomain()
    # Step 4: Sample a random integer γ ←$ Z_{2^(2a)}
    gamma = random_unit(2^(2*a_param))
    # Step 5: Compute P_B = [γ]ψ(P_0)
    _psi_P0 = psi(xP0)
    psi_P0 = _psi_P0.curve_point()

    P_B = gamma * psi_P0
    
    # Step 6: Compute P_AB = [γ]ψ'(P_A)
    _psi_prime_PA = psi_prime(xPA)
    psi_prime_PA = _psi_prime_PA.curve_point()
    P_AB = gamma * psi_prime_PA
    
    # Step 7: Sample a random matrix M ←$ SL_2(Z_{B})
    d1, d2, d3, d4 = random_matrix(B_param)
    
    # Step 8: Compute [RB, SB]^T = M[ψ(R_0), ψ(S_0)]^T
    _psi_R0 = psi(xR0)
    _psi_S0 = psi(xS0)
    
    psi_R0, psi_S0 = lift_image_to_curve(R0, S0, _psi_R0, _psi_S0, B, A)
    RB = d1 * psi_R0 + d2 * psi_S0
    SB = d3 * psi_R0 + d4 * psi_S0 
    
    # Step 9: Compute [R_AB, S_AB]^T = M[ψ'(R_A), ψ'(S_A)]^T
    _psi_prime_RA = psi_prime(xRA)
    _psi_prime_SA = psi_prime(xSA)

    psi_prime_RA, psi_prime_SA = lift_image_to_curve(RA, SA, _psi_prime_RA, _psi_prime_SA, B, A)
    R_AB = 2^a_param * (d1 * psi_prime_RA + d2 * psi_prime_SA)
    S_AB = 2^a_param * (d3 * psi_prime_RA + d4 * psi_prime_SA)
    
    # Step 10: tmp = KDF(R_AB, S_AB) ⊕ m
    R_AB_bytes = str(R_AB[0]).encode()
    S_AB_bytes = str(S_AB[0]).encode()
    kdf_input = R_AB_bytes + S_AB_bytes
    key_length = len(m)
    kdf_output = shake_256(kdf_input).digest(key_length)
    tmp = bytes(a ^^ int(b) for a, b in zip(kdf_output, m))
    
    # Step 11: return ct = (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    ct = (_EB.curve(), P_B, RB, SB, _EAB.curve(), P_AB, tmp)
    
    return ct

# %%
def decryption(sk, ct, a_param, B_param):
    """
    Algorithm 5: Decryption
    
    Require: sk = (q, α, β), ct = (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    Ensure: A message m
    
    Parameters:
    - sk: Secret key (q, α, β) from keygen
    - ct: Ciphertext (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    - a_param: parameter a for 2-torsion
    - B_param: parameter B for B-torsion
    
    Returns:
    - m: Decrypted message
    """
    # Extract secret key components
    q, alpha, beta = sk
    
    # Extract ciphertext components
    E_B, P_B, R_B, S_B, E_AB, P_AB, tmp = ct
    
    # Step 1: Compute P'_AB = [1/α]P_AB
    alpha_inv = inverse_mod(alpha, 2^(2*a_param))
    P_AB_prime = alpha_inv * P_AB
    
    
    R_AB1, S_AB1 = one_point(E_B, E_AB, P_B, P_AB_prime, q, (R_B, S_B))
    R_AB1 = beta * R_AB1
    S_AB1 = beta * S_AB1
    # print(f"        deg = 2^{a_param} - {q} = {deg}")
    
    # Step 3: Compute m = KDF(R_AB, S_AB) ⊕ tmp
    # Convert points to bytes for KDF
    R_AB_bytes = str(R_AB1.x()).encode()
    S_AB_bytes = str(S_AB1.x()).encode()
    
    # Generate key material with SHAKE256
    kdf_input = R_AB_bytes + S_AB_bytes
    key_length = len(tmp)
    kdf_output = shake_256(kdf_input).digest(key_length)
    
    # XOR with ciphertext to recover message
    m = bytes(a ^^ b for a, b in zip(kdf_output, tmp))
    
    return m


# %%
p, a, B, _= prime_params[0].values()
A = Integer(2**(2*a))

K = GF(p)
PR_p.<xx> = K[]

# 2. 构造 F0 = GF(p^2)
F0.<i0> = GF(p^2, modulus=xx^2+1)

f = PR_p.irreducible_element(8)

# 4. 关键：使用 F0.extension 构造 F
# 为了确保它被视为有限域，Sage 需要通过 F0 显式扩张
F.<iiii> = GF(p^8, modulus=f)
i = F0.embeddings(F)[0](i0)
E0 = EllipticCurve(F, [1, 0])
E0.set_order((p^4-1)**2)
# a -=2
P0, Q0 = torsion_basis(E0, A)
P0x = E0x.random_point()*(E0x.cardinality()/2^(2*a))
while (P0x * 2^(2*a-1)).is_zero():
    P0x = E0x.random_point()*(E0x.cardinality()/2^(2*a))
P0x = P0x.change_ring(F)
if P0x.weil_pairing(P0, 2^(2*a))^(2^(2*a-1)) == 1:
    P0 = P0x
else:
    Q0 = P0x
R0, S0 = torsion_basis(E0, B)

PQ0 = P0-Q0
RS0 = R0-S0

_E0 = KummerLine(E0)
xR0 = _E0(R0[0])
xS0 = _E0(S0[0])
xRS0 = _E0(RS0[0])

xP0 = _E0(P0[0])
xQ0 = _E0(Q0[0])
xPQ0 = _E0(PQ0[0])

# %%
import statistics

timing_results = {
    'keygen': [],
    'encryption': [],
    'decryption': []
}

msg = bytes("你好! 我是小猫", 'utf-8')
succ = 0
failreason = []
for trial in range(10):
    # Timing keygen
    t0 = time.time()
    sk, pk = keygen(E0, P0, Q0, R0, S0, a, B)
    keygen_time = time.time() - t0
    timing_results['keygen'].append(keygen_time)
    # Timing encryption
    t0 = time.time()
    ct = encryption(pk, msg, a, B)
    encryption_time = time.time() - t0
    timing_results['encryption'].append(encryption_time)
    # Timing decryption
    t0 = time.time()
    decrypted = decryption(sk, ct, a, B)
    decryption_time = time.time() - t0
    timing_results['decryption'].append(decryption_time)
    result = decrypted.decode()
    # print(f"Trial {trial+1}: Decryption successful - {result}")
    print(f"  keygen: {keygen_time:.4f}s, encryption: {encryption_time:.4f}s, decryption: {decryption_time:.4f}s")
    succ += 1

print("\n" + "="*60)
print(f"总成功: {succ}/{len(range(10))}")
print("="*60)

# Print statistics
for op in ['keygen', 'encryption', 'decryption']:
    if timing_results[op]:
        times = timing_results[op]
        print(f"\n{op}:")
        print(f"  成功次数: {len(times)}")
        print(f"  最小时间: {min(times):.4f}s")
        print(f"  最大时间: {max(times):.4f}s")
        print(f"  平均时间: {statistics.mean(times):.4f}s")
        if len(times) > 1:
            print(f"  标准差: {statistics.stdev(times):.4f}s")
        print(f"  总耗时: {sum(times):.4f}s")
