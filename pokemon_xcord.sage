# %%
import os, gc
import random
import time
from hashlib import shake_256
from sage.all import *
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from montgomery_isogenies.isogenies_x_only import lift_image_to_curve

from utilities.discrete_log import discrete_log_pari
from utilities.supersingular import torsion_basis, torsion_basis_2e
from utilities.utils import speed_up_sagemath
from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_endomorphism import KaniEndoHalf
from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
gc.disable()
speed_up_sagemath()
prime_params = [{'prime': 38834938483610889734772360524892636617924059316004375193980790404642840335381409451474943,
 'c': 19,
 'a': 126,
 'b': 24,
 'bits': 295}, {'prime': 1437109229254356717681521168159299154183356802192981571221945465128594431010067536245211284485387759995646851969299095099429251971022847,
 'c': 12,
 'a': 191,
 'b': 40,
 'bits': 449}, {'prime': 300794947883270440918790320272237309765137040234666071462510015573632691421653784724779362035035310704193643738753950239846195422171320926207701270866078070654007289891491151871,
 'c': 18,
 'a': 253,
 'b': 48,
 'bits': 587}]

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
    F = factor(n) 
    for p, exponent in F:
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
def search_special_prime(lam, c_min=1, c_max=20, iterations=100):
    
    # a ~ lam
    # b * log(3) ~ (lam / 3) * log(2)  => b ~ (lam * log(2)) / (3 * log(3))
    
    log2 = float(log(2))
    log3 = float(log(3))
    
    target_b = int((lam * log2) / (3 * log3))
    
    a_min = lam - 1
    a_max = lam + 1
    
    b_min = max(1, target_b - 5)
    b_max = target_b + 5

    for _ in range(iterations):
        curr_c = randint(c_min, c_max)
        curr_a = randint(a_min, a_max)
        curr_b = randint(b_min, b_max)
        
        # n = c * 2^(2*a) * 3^b - 1
        candidate = curr_c * (2**(2 * curr_a)) * (3**curr_b) - 1
        
        if candidate % 2 == 0:
            continue
        
        if candidate.is_pseudoprime():
            if candidate.is_prime():
                return {
                    "prime": candidate,
                    "c": curr_c,
                    "a": curr_a,
                    "b": curr_b,
                    "bits": candidate.nbits()
                }
            
    return None

# %%
def gen_isogeny_deg(E0, P0, Q0,R0,S0, q, a_param, b_param):
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
    - b_param: parameter b for 3-torsion
    
    Returns:
    - E1: Codomain curve
    - P1, Q1: Torsion basis of E_1[2^(2a)]
    """
    # Step 1: Generate an endomorphism θ of degree q(2^(2a) - q)3^b
    deg_theta = q * (2^(2*a_param) - q) * (3^b_param)

    
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
        wp_ = R0_prime.weil_pairing(R, 3^b_param, algorithm='pari')
        wp = S0_prime.weil_pairing(R, 3^b_param, algorithm='pari')
        
        discrete_log_pari(wp_, wp, 3^b_param)

        K3_dual = S0_prime
    except TypeError:
        R = R0_prime
        K3_dual = R0_prime # P0_prime is guaranteed to be lin. indep with K3
        
    phi = E0.isogeny(K3_dual, algorithm="factored", model="montgomery")
    
    
    # Step 5: Compute P_1, Q_1 = [1/3^b]φ(P'_0), [1/3^b]φ(Q'_0)
    # Evaluate phi on the endomorphism images
    phi_P0_prime = phi(P0_prime)
    phi_Q0_prime = phi(Q0_prime)
    
    # Divide by 3^b (multiply by the inverse)
    inv_3b = pow(3^b_param, -1, 2^(2*a_param))  # Inverse in the group order
    P1 = inv_3b * phi_P0_prime
    Q1 = inv_3b * phi_Q0_prime
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
def one_point(E0, E1, P0_op, phi1P0, d, points, tempx, tempy):
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
    # a1, a2 = two_squares(n-d)
    a1, a2 = tempx, tempy
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
    F = KaniEndoHalf(P2, Q2, P3, (2*int(P2.weil_pairing(Q2, 2^a, algorithm='pari')^d == P3.weil_pairing(Q3, 2^a,algorithm='pari'))-1)*Q3, d, a1, a2, a, a)
    phi1prime = lambda pt: F(TuplePoint(pt, E2(0), E3(0), E3(0)))[2]
    return map(lambda pt: phi2_prime_dual(xE3(phi1prime(phi2(_E0(pt)).curve_point()))), points)

# %%
def keygen(E0, P0, Q0, R0, S0, a_param, b_param):
    """
    Algorithm 3: KeyGen
    
    Require: E_0, E_0[2^(2a)] = ⟨P_0, Q_0⟩, E_0[3^b] = ⟨R_0, S_0⟩.
    Ensure: sk, pk.
    
    Parameters:
    - E0: Supersingular elliptic curve
    - P0, Q0: Basis of E_0[2^(2a)]
    - R0, S0: Basis of E_0[3^b]
    - a_param: parameter a for 2-torsion
    - b_param: parameter b for 3-torsion
    
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
    
    phi = gen_isogeny_deg(E0,P0,Q0,R0,S0,q, a_param, b_param)
    
    P23x = P0+R0
    Q23x = Q0+S0
    alpha = random_unit(2^(2*a_param))
    beta = random_unit(3^b_param)
    imP23x = phi(P23x)
    imQ23x = phi(Q23x)
    imPQ23x = phi(P23x+Q23x)
    if not (imP23x+imQ23x == imPQ23x or imP23x+imQ23x == -imPQ23x):
        imQ23x = -imQ23x
    EA = imP23x.curve()
    _EA = KummerLine(EA)
    PA = alpha * inverse_mod(3^b_param, two_to_twoa) * (3^b_param * imP23x)
    QA = alpha * inverse_mod(3^b_param, two_to_twoa) * (3^b_param * imQ23x)
    
    RA = beta * inverse_mod(two_to_twoa, 3^b_param) * (two_to_twoa * imP23x)
    SA = beta * inverse_mod(two_to_twoa, 3^b_param) * (two_to_twoa * imQ23x)
    
    _, xPA, xQA, xPQA = point_to_xonly(PA, QA)
    _, xRA, xSA, xRSA = point_to_xonly(RA, SA)

    sk = (q, alpha, beta, tempx, tempy)
    pk = (_EA, xPA, xQA, xPQA, xRA, xSA, xRSA, RA, SA)
    
    
    return sk, pk

# %%
def encryption(pk, m, a_param, b_param):
    """
    Algorithm 4: Encryption
    
    Require: A message m, a public key pk = (E_A, xPA, xQA, xPQA, xRA, xSA, xRSA).
    Ensure: A ciphertext ct.
    
    Parameters:
    - pk: Public key (E_A, xPA, xQA, xPQA, xRA, xSA, xRSA) from keygen
    - m: Message to encrypt (bytes)
    - a_param: parameter a for 2-torsion
    - b_param: parameter b for 3-torsion
    
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
    
    # Step 7: Sample a random matrix M ←$ SL_2(Z_{3^b})
    d1, d2, d3, d4 = random_matrix(3^b_param)
    
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
def decryption(sk, ct, a_param, b_param):
    """
    Algorithm 5: Decryption
    
    Require: sk = (q, α, β), ct = (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    Ensure: A message m
    
    Parameters:
    - sk: Secret key (q, α, β) from keygen
    - ct: Ciphertext (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    - a_param: parameter a for 2-torsion
    - b_param: parameter b for 3-torsion
    
    Returns:
    - m: Decrypted message
    """
    # Extract secret key components
    q, alpha, beta, tempx, tempy = sk
    
    # Extract ciphertext components
    E_B, P_B, R_B, S_B, E_AB, P_AB, tmp = ct
    
    # Step 1: Compute P'_AB = [1/α]P_AB
    alpha_inv = inverse_mod(alpha, 2^(2*a_param))
    P_AB_prime = alpha_inv * P_AB
    
    
    R_AB1, S_AB1 = one_point(E_B, E_AB, P_B, P_AB_prime, q, (R_B, S_B), tempx, tempy)
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
for pval in prime_params:
    p, c, a, b, _= pval.values()
    A = Integer(2**(2*a))
    B = Integer(3**b)
    FFF = GF(p)
    FF.<xx> = GF(p)[]
    F.<i> = GF(p^2, modulus=xx^2+1)
    E0 = EllipticCurve(F, [1, 0])
    E0.set_order((p+1)**2)
    E0x = E0.change_ring(FFF)
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
    run_time = 100

    for trial in range(run_time):
        # Timing keygen
        t0 = time.time()
        sk, pk = keygen(E0, P0, Q0, R0, S0, a, b)
        keygen_time = time.time() - t0
        timing_results['keygen'].append(keygen_time)
        # Timing encryption
        t0 = time.time()
        ct = encryption(pk, msg, a, b)
        encryption_time = time.time() - t0
        timing_results['encryption'].append(encryption_time)
        # Timing decryption
        t0 = time.time()
        decrypted = decryption(sk, ct, a, b)
        decryption_time = time.time() - t0
        timing_results['decryption'].append(decryption_time)
        result = decrypted.decode()
        # print(f"Trial {trial+1}: Decryption successful - {result}")
        # print(f"  keygen: {keygen_time:.4f}s, encryption: {encryption_time:.4f}s, decryption: {decryption_time:.4f}s")
        succ += 1

    print("\n" + "="*60)
    print(f"{p.nbits()}, success: {succ}/{run_time}")
    print("="*60)

    # Print statistics
    for op in ['keygen', 'encryption', 'decryption']:
        if timing_results[op]:
            times = timing_results[op]
            print(f"\n{op}:")
            print(f"  success count: {len(times)}")
            print(f"  minimum: {min(times):.4f}s")
            print(f"  maximum: {max(times):.4f}s")
            print(f"  average: {statistics.mean(times):.4f}s")
            if len(times) > 1:
                print(f"  stddev: {statistics.stdev(times):.4f}s")
