# %%
import os
import random
import time
from hashlib import shake_256
from IPython import embed
from sage.all import pari
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from montgomery_isogenies.isogenies_x_only import lift_image_to_curve

from utilities.discrete_log import BiDLP, discrete_log_pari
from utilities.supersingular import torsion_basis, torsion_basis_with_pairing

from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_endomorphism import KaniEndoHalf
from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
from sage.rings.sum_of_squares import is_sum_of_two_squares_pyx
proof.all(False)

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
# %%

def manual_factored_isogeny(E, P, factors):
    """
    E: 起始曲线
    P: kernel 生成点
    factors: 列表，包含所有素因子 (含重复)，例如 [2, 2, 2, 3, 3]
    """
    curr_E = E
    curr_P = P
    phi_total = None
    
    # 拷贝一份因子列表，避免修改原数据
    remaining_factors = list(factors)
    
    while remaining_factors:
        # 1. 提取当前的素因子 l
        l = remaining_factors.pop(0)
        
        # 2. 计算余项乘积 M = \prod \ell_j
        M = 1
        for f in remaining_factors:
            M *= f
            
        # 3. 找到阶为 l 的点 Q
        Q = M * curr_P
        
        if Q.is_zero():
            # 如果 Q 是 0，说明 P 的阶中不含当前的 l，直接跳过
            continue
            
        # 4. 构造小步长同源
        # 明确传入 degree=l，强制 Sage 使用 Velu 而不触发 order() 检查
        phi = curr_E.isogeny(Q, degree=l, check=False)
        
        # 5. 更新状态
        curr_E = phi.codomain()
        curr_P = phi(curr_P)
        
        if phi_total is None:
            phi_total = phi
        else:
            phi_total = phi * phi_total
            
    return phi_total
# 给定一个偏移量（seed），寻找第一个合法的点
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

# 示例使用
# P_basis, Q_basis = get_torsion_basis_fast(E8, 2^3)
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
    - B_param: parameter B for B-torsion
    
    Returns:
    - E1: Codomain curve
    - P1, Q1: Torsion basis of E_1[2^(2a)]
    """
    # Step 1: Generate an endomorphism θ of degree q(2^(2a) - q)B
    deg_theta = q * (2^(2*a_param) - q) * B_param

    
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
        print(sq)
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
    theta_dual = lambda P: xx*P - yy*i_end(P) - zz*pi_end(P) - tt*i_end(pi_end(P))
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
        wp_ = R0_prime.weil_pairing(R, B_param)
        wp = S0_prime.weil_pairing(R, B_param)
        
        discrete_log_pari(wp_, wp, B_param)

        K3_dual = S0_prime
    except TypeError:
        R = R0_prime
        K3_dual = R0_prime # P0_prime is guaranteed to be lin. indep with K3
    
    # phi = E0.isogeny(K3_dual, algorithm="factored", model="montgomery")
    phi = manual_factored_isogeny(E0, K3_dual, [p for p, e in B_param.factor() for _ in range(e)])
    print("phi generated")
    
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
        e_P0_Q0 = P0.weil_pairing(Q0, 2^(2*a_param))
        e_P2_Q2 = P2.weil_pairing(Q2, 2^(2*a_param))
        e_P3_Q3 = P3.weil_pairing(Q3, 2^(2*a_param))
        
        target_pairing = e_P0_Q0^q
        if e_P2_Q2 == target_pairing:
            print(f"Step 9: Pairing check passed, returning (E_2, P_2, Q_2)")
            return prodeval(Phi, E1(0), 0)
        else:
            assert e_P3_Q3 == target_pairing
            print(f"Step 9: Pairing check failed, returning (E_3, P_3, Q_3)")
            return prodeval(Phi, E1(0), 1)
    except:
      raise RuntimeError("Weil pairing computation failed.")

# %%
def one_point(E0, E1, P0, phi1P0, d, points):
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
    
    P1 = phi1P0
    n = 2^a  # degree for auxiliary isogenies
    a1, a2 = two_squares(n-d)
    if P0.order()<2^(2*a) or P1.order()<2^(2*a):
        raise Exception("P0, P1 not order 2^2a")
    Q0_basis = get_torsion_basis_fast(E0, 2^(2*a), (p^4-1), random=False)
    # Try to find Q0 such that the pairing is non-trivial
    # This should be changed to a random search in practice
    zeta = None
    for i in range(0, 1000):
        Q0_cand = Q0_basis[0] + i * Q0_basis[1]  # Try multiples of the second basis point
        try:
            pairing_val = P0.weil_pairing(Q0_cand, 2^(2*a))
            if pairing_val^(2^(2*a-1)) != 1:
                Q0 = Q0_cand
                zeta = pairing_val
                break
        except:
            continue
    
    if zeta is None:
        raise ValueError("Could not find suitable Q0")
    K2 = 2^a * P0  # Kernel generator
    
    # Use isogeny method directly
    phi2 = E0.isogeny(K2, algorithm="factored", model="montgomery")
    E2 = phi2.codomain()
    # Compute torsion points (P2, Q2) = (φ2(P0), 2^a φ2(Q0))
    P2 = phi2(P0)
    Q0_image = phi2(Q0)
    Q2 = 2^a * Q0_image  # 2^a φ2(Q0)
    
    # Step 3: Select Q1 ∈ E1[2^(2a)] such that e_{2^(2a)}(P1, Q1) = ζ^d
    Q1_basis = get_torsion_basis_fast(E1, 2^(2*a), (p^4-1), random=False)
    for i in range(0, 2^(2*a)):
        Q1_cand = Q1_basis[0] + i * Q1_basis[1]  # Try multiples of the second basis point
        try:
            pairing_val = P1.weil_pairing(Q1_cand, 2^(2*a))
            if pairing_val^(2^(2*a-1)) != 1:
                S1 = Q1_cand
                break
        except:
            continue
    target_pairing = zeta^d
    ee = P1.weil_pairing(S1, 2^(2*a))
    ii = discrete_log_pari(target_pairing, ee, 2^(2*a))
    Q1 = P1+ii*S1
    
    K2_prime = 2^a * P1
    
    # Use isogeny method directly
    phi2_prime = E1.isogeny(K2_prime, algorithm="factored", model="montgomery")
    E3 = phi2_prime.codomain()
    P3 = phi2_prime(P1)
    Q1_image = phi2_prime(Q1)
    Q3 = 2^a * Q1_image
    
    # 原有的代码继续...
    F = KaniEndoHalf(P2, Q2, P3, Q3, d, a1, a2, a, a)
    phi1prime = lambda pt: F(TuplePoint(pt, E2(0), E3(0), E3(0)))[2]
    return map(lambda pt: phi2_prime.dual()(phi1prime(phi2(pt))), points)

# %%
def keygen(E0, P0, Q0, R0, S0, a_param, B_param):
    """
    Algorithm 3: KeyGen
    
    Require: E_0, E_0[2^(2a)] = ⟨P_0, Q_0⟩, E_0[3^b] = ⟨R_0, S_0⟩.
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
        tempx = ZZ(random.randint(1, two_to_a.isqrt())//2)*2
        tempy = ZZ(random.randint(1, two_to_a.isqrt())//2)*2+1
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
    
    PA = alpha * inverse_mod(B_param, two_to_twoa) * (B_param * imP23x)
    QA = alpha * inverse_mod(B_param, two_to_twoa) * (B_param * imQ23x)
    
    RA = beta * inverse_mod(two_to_twoa, B_param) * (two_to_twoa * imP23x)
    SA = beta * inverse_mod(two_to_twoa, B_param) * (two_to_twoa * imQ23x)
    
    sk = (q, alpha, beta)
    pk = (EA, PA, QA, RA, SA)
    
    
    return sk, pk

# %%
def encryption(pk, m, a_param, B_param):
    """
    Algorithm 4: Encryption
    
    Require: A message m, a public key pk = (E_A, P_A, Q_A, R_A, S_A).
    Ensure: A ciphertext ct.
    
    Parameters:
    - pk: Public key (E_A, P_A, Q_A, R_A, S_A) from keygen
    - m: Message to encrypt (bytes)
    - a_param: parameter a for 2-torsion
    - B_param: parameter B for B-torsion
    
    Returns:
    - ct: Ciphertext (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    """
    E_A, P_A, Q_A, R_A, S_A = pk
    
    # Step 1: Select a random integer s ∈ [0, 2^(2a)]
    s = ZZ(random.randint(0, 2^(2*a_param-1)))*2+1
    kernel_B = P0 + s * Q0
    # Step 2: Compute the isogeny ψ : E_0 → E_B with kernel ⟨P_0 + sQ_0⟩
    
    while (2^(2*a_param-1)*kernel_B)[0] == 0:
        s = ZZ(random.randint(0, 2^(2*a_param-1)))*2+1
        kernel_B = P0 + s * Q0

    psi = E0.isogeny(kernel_B, algorithm="factored", model="montgomery")
    
    E_B = psi.codomain()
    
    # Step 3: Compute the isogeny ψ' : E_A → E_AB with kernel ⟨P_A + sQ_A⟩
    kernel_AB = P_A + s * Q_A
    psi_prime = E_A.isogeny(kernel_AB, algorithm="factored", model="montgomery")
    
    E_AB = psi_prime.codomain()
    
    # Step 4: Sample a random integer γ ←$ Z_{2^(2a)}
    gamma = random_unit(2^(2*a_param))
    # Step 5: Compute P_B = [γ]ψ(P_0)
    psi_P0 = psi(P0)
    P_B = gamma * psi_P0
    
    # Step 6: Compute P_AB = [γ]ψ'(P_A)
    psi_prime_PA = psi_prime(P_A)
    P_AB = gamma * psi_prime_PA
    
    # Step 7: Sample a random matrix M ←$ SL_2(Z_{B})
    d1, d2, d3, d4 = random_matrix(B_param)
    
    # Step 8: Compute [R_B, S_B]^T = M[ψ(R_0), ψ(S_0)]^T
    psi_R0 = psi(R0)
    psi_S0 = psi(S0)
    
    R_B = d1 * psi_R0 + d2 * psi_S0
    S_B = d3 * psi_R0 + d4 * psi_S0 
    
    # Step 9: Compute [R_AB, S_AB]^T = M[ψ'(R_A), ψ'(S_A)]^T
    psi_prime_RA = psi_prime(R_A)
    psi_prime_SA = psi_prime(S_A)

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
    ct = (E_B, P_B, R_B, S_B, E_AB, P_AB, tmp)
    
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
    R_AB_bytes = str(R_AB1[0]).encode()
    S_AB_bytes = str(S_AB1[0]).encode()
    
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
print("get P0 Q0")
P0, Q0 = get_torsion_basis_fast(E0, A, (p^4-1), random=False)
# P0, Q0 = E0.torsion_basis(A)
print("get R0 S0")
R0, S0 = get_torsion_basis_fast(E0, B, (p^4-1), random=False)
# R0, S0 = E0.torsion_basis(B)



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
    # try:
        # Timing keygen
        t0 = time.time()
        sk, pk = keygen(E0, P0, Q0, R0, S0, a, B)
        keygen_time = time.time() - t0
        timing_results['keygen'].append(keygen_time)
        print("Keygen completed.")
        # Timing encryption
        t0 = time.time()
        ct = encryption(pk, msg, a, B)
        encryption_time = time.time() - t0
        timing_results['encryption'].append(encryption_time)
        print("Encryption completed.")
        # Timing decryption
        try:
            t0 = time.time()
            decrypted = decryption(sk, ct, a, B)
            decryption_time = time.time() - t0
            timing_results['decryption'].append(decryption_time)
            
            result = decrypted.decode()
            print(f"Trial {trial+1}: Decryption successful - {result}")
            print(f"  keygen: {keygen_time:.4f}s, encryption: {encryption_time:.4f}s, decryption: {decryption_time:.4f}s")
            succ += 1
        except Exception as e:
            print(f"Trial {trial+1}: Decryption failed - {e}")
            
    # except Exception as e:
    #     failreason.append(str(e))
    #     print(f"Trial {trial+1}: Error - {e}")

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
