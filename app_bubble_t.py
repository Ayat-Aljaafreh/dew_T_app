import streamlit as st
import math

st.title("Bubble T Calculation")

# ===================== Utility Functions =====================

def calculate_T0(P, A, B, C, P_unit_Ant, T_unit_Ant,x):
    Ti_sat = []                                     #Calculate Tsat
    n = len(A)

    for i in range(n):
        T_i = B[i] / (A[i] - math.log(P)) - C[i]  # in same unit as Antoine
        Ti_sat.append(T_i)

    T0 = 0.0
    for i in range(n):
        T0 += x[i] * Ti_sat[i]

    return T0


def calculate_Psat(T_1, T_unit_Ant, A, B, C, P_unit_Ant):
    Psat = []
    n = len(A)

    for i in range(n):

        Psat_i = math.exp((A[i] - B[i] / (T_1 + C[i])))
        Psat.append(Psat_i)

    return Psat



def calculate_Lambda(a, V, R, T_Ant, T_unit):

    T_K = to_kelvin(T_Ant, T_unit)
    n = len(V)
    Lambda = [[0.0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                Lambda[i][j] = 1.0
            else:
                Lambda[i][j] = (V[j] / V[i]) * math.exp(-a[i][j] / (R * T_K))

    return Lambda



def calculate_gamma(x, Lambda):
    n = len(x)
    gamma = []

    for i in range(n):

        # first sum: sum_j (xj * Lambda_ij)
        sum1 = 0.0
        for j in range(n):
            sum1 += x[j] * Lambda[i][j]

        # second term: sum_k [ xk * Lambda_ki / sum_j(xj * Lambda_kj) ]
        sum2 = 0.0
        for k in range(n):
            denom = 0.0
            for j in range(n):
                denom += x[j] * Lambda[k][j]

            sum2 += x[k] * Lambda[k][i] / denom

        ln_gamma_i = 1.0 - math.log(sum1) - sum2
        gamma.append(math.exp(ln_gamma_i))

    return gamma


def convert_pressure_to_kpa(P, unit):
    if unit == "kpa":
        return P
    elif unit == "bar":
        return P * 100.0
    elif unit == "mmhg":
        return P * 0.133322
    else:
        raise ValueError("Unsupported pressure unit")

def to_kelvin(T, unit):
    if unit == "k":
        return T
    elif unit == "c":
        return T + 273
    else:
        raise ValueError("Invalid temperature unit")

def calculate_phi(T_old, P, Bik, R, yi, T_unit):

    T_K = to_kelvin(T_old, T_unit)

    n = len(yi)
    phi = []

    # Compute Sij = 2*Bik - Bii - Bkk
    S = [[2 * Bik[i][j] - Bik[i][i] - Bik[j][j] for j in range(n)] for i in range(n)]

    for k in range(n):
        sum_term = 0.0
        for i in range(n):
            for j in range(n):
                sum_term += yi[i] * yi[j] * (2 * S[i][k] - S[i][j])

        ln_phi_k = (P / (R * T_K)) * (Bik[k][k] + 0.5 * sum_term)
        phi_k = math.exp(ln_phi_k)
        phi.append(phi_k)

    return phi





def calculate_alpha(Pi_sat, k_index):
    alpha = [Pi_sat[i] / Pi_sat[k_index] for i in range(len(Pi_sat))]
    return alpha


# ===================== Bubble Point Iteration =====================

def bubble_point_iteration(T0, P, x, A, B, C, T_unit_Ant, V, a, R,
                           k_index, P_unit_Ant, tol,
                           nonideal_liq, nonideal_gas, Bik, max_iter=100):

    T_old = T0
    iteration = 0
    phi = [1.0] * len(x)

    Pi_sat = calculate_Psat(T_old, T_unit_Ant, A, B, C, P_unit_Ant)

    if nonideal_liq:
        Lambda = calculate_Lambda(a, V, R, T_old, T_unit_Ant)
        gamma = calculate_gamma(x, Lambda)
    else:
        gamma = [1.0] * len(x)

    alpha = calculate_alpha(Pi_sat, k_index)
    numerator_sum = sum(((x[i] * gamma[i] * alpha[i]) / phi[i]) for i in range(len(x)))
    Pk_sat = P / numerator_sum
    Tk_sat = B[k_index] / (A[k_index] - math.log(Pk_sat)) - C[k_index]

    while iteration < max_iter:

        iteration += 1
        T_old = Tk_sat

        Pi_sat = calculate_Psat(T_old, T_unit_Ant, A, B, C, P_unit_Ant)

        # ---- Liquid phase ----
        if nonideal_liq:
            Lambda = calculate_Lambda(a, V, R, T_old, T_unit_Ant)
            gamma = calculate_gamma(x, Lambda)
        else:
            gamma = [1.0] * len(x)

        yi = [(x[i] * gamma[i] * Pi_sat[i]) / (P * phi[i]) for i in range(len(x))]

        # ---- Gas phase ----
        if nonideal_gas:
            phi = calculate_phi(T_old, P, Bik, R, yi, T_unit_Ant)
        else:
            phi = [1.0] * len(x)

        alpha = calculate_alpha(Pi_sat, k_index)
        numerator_sum = sum(((x[i] * gamma[i] * alpha[i]) / phi[i]) for i in range(len(x)))
        Pk_sat = P / numerator_sum
        Tk_sat = B[k_index] / (A[k_index] - math.log(Pk_sat)) - C[k_index]

        error = (abs(T_old - Tk_sat) / T_old) * 100

        # ===== Streamlit Output =====
        st.text(f"\nIteration {iteration}")
        st.text(f"T_old = {T_old:.8f} {T_unit_Ant.upper()}")
        st.text(f"Tk_sat = {Tk_sat:.8f} {T_unit_Ant.upper()}")
        st.text(f"Error  = {error:.8f} %")
        st.text(f"Pk_sat = {Pk_sat:.8f} kPa")

        st.text("\nComponent-wise properties:")
        st.text("i   Pi_sat (kPa)     gamma_i        alpha_ik")
        st.text("-" * 60)

        for i in range(len(x)):
            st.text(f"{i+1:<3} {Pi_sat[i]:<14.8f} {gamma[i]:<14.8f} {alpha[i]:<14.8f}")

        if error <= tol:
            T_final = Tk_sat
            break

    else:
        T_final = Tk_sat

    Pi_sat = calculate_Psat(T_final, T_unit_Ant, A, B, C, P_unit_Ant)

    if nonideal_liq:
        Lambda = calculate_Lambda(a, V, R, T_final, T_unit_Ant)
        gamma = calculate_gamma(x, Lambda)
    else:
        gamma = [1.0] * len(x)

    yi = [(x[i] * gamma[i] * Pi_sat[i]) / (P * phi[i]) for i in range(len(x))]

    st.success("\nConverged Bubble Point Results:")
    st.text(f"Final T = {T_final:.8f} {T_unit_Ant.upper()}, Iterations = {iteration}")
    st.text("Vapor phase mole fractions yi:")

    for i, y in enumerate(yi):
        st.text(f"y{i+1} = {y:.8f}")

    return T_final, yi


# ===================== Streamlit UI =====================

# ===================== Streamlit UI =====================

n = st.sidebar.number_input("Number of components", min_value=1, value=2, key="n")

nonideal_liq = st.sidebar.radio("Liquid phase non-ideal?", ["Yes", "No"], index=0, key="nonideal_liq") == "Yes"
nonideal_gas = st.sidebar.radio("Gas phase non-ideal?", ["Yes", "No"], index=0, key="nonideal_gas") == "Yes"

# Wilson parameters
a = [[0.0]*n for _ in range(n)]
V = []

if nonideal_liq:
    st.sidebar.subheader("Wilson Parameters")
    for i in range(n):
        for j in range(n):
            if i != j:
                key_a = f"a{i}{j}"
                if key_a not in st.session_state:
                    st.session_state[key_a] = 0.0
                a[i][j] = st.sidebar.number_input(f"a{i+1}{j+1}", value=st.session_state[key_a], format="%.8f", key=key_a)

    for i in range(n):
        key_V = f"V{i}"
        if key_V not in st.session_state:
            st.session_state[key_V] = 1.0
        V_i = st.sidebar.number_input(f"V{i+1}", value=st.session_state[key_V], format="%.8f", key=key_V)
        V.append(V_i)
else:
    V = [1.0]*n

# Virial coefficients
Bik = [[0.0]*n for _ in range(n)]
if nonideal_gas:
    st.sidebar.subheader("Virial Coefficients")
    for i in range(n):
        for j in range(i, n):
            key_Bik = f"Bik{i}{j}"
            if key_Bik not in st.session_state:
                st.session_state[key_Bik] = 0.0
            val = st.sidebar.number_input(f"B{i+1}{j+1}", value=st.session_state[key_Bik], format="%.8f", key=key_Bik)
            Bik[i][j] = Bik[j][i] = val

# Antoine constants
A = []
B = []
C = []
for i in range(n):
    keyA = f"A{i}"
    keyB = f"B{i}"
    keyC = f"C{i}"
    
    if keyA not in st.session_state:
        st.session_state[keyA] = 8.07131
    if keyB not in st.session_state:
        st.session_state[keyB] = 1730.63
    if keyC not in st.session_state:
        st.session_state[keyC] = 233.426

    Ai = st.sidebar.number_input(f"A{i+1}", value=st.session_state[keyA], format="%.8f", key=keyA)
    Bi = st.sidebar.number_input(f"B{i+1}", value=st.session_state[keyB], format="%.8f", key=keyB)
    Ci = st.sidebar.number_input(f"C{i+1}", value=st.session_state[keyC], format="%.8f", key=keyC)

    A.append(Ai)
    B.append(Bi)
    C.append(Ci)

# System parameters
if "Ptot" not in st.session_state:
    st.session_state.Ptot = 101.325
Ptot = st.sidebar.number_input("System Pressure", value=st.session_state.Ptot, format="%.8f", key="Ptot")

if "R" not in st.session_state:
    st.session_state.R = 8.314
R = st.sidebar.number_input("Gas Constant R", value=st.session_state.R, format="%.8f", key="R")

# Mole fractions
x = []
for i in range(n):
    key_x = f"x{i}"
    if key_x not in st.session_state:
        st.session_state[key_x] = 1/n
    xi = st.sidebar.number_input(f"y{i+1}", value=st.session_state[key_x], format="%.8f", key=key_x)
    x.append(xi)

if "tol" not in st.session_state:
    st.session_state.tol = 0.01
tol = st.sidebar.number_input("Temperature tolerance (%)", value=st.session_state.tol, format="%.8f", key="tol")

# Temperature and Pressure units
T_unit = st.sidebar.selectbox("Temperature unit", ["C", "K"], index=0, key="T_unit").lower()
P_unit = st.sidebar.selectbox("Pressure unit", ["mmHg", "bar", "kPa"], index=2, key="P_unit").lower()


# ===================== Run =====================

if st.button("Calculate Bubble Point"):
    T0 = calculate_T0(P, A, B, C, P_unit_Ant, T_unit_Ant, x)

    T_final, y_final = bubble_point_iteration(
        T0, P, x, A, B, C, T_unit_Ant,
        V, a, R, k_index,
        P_unit_Ant, tol,
        nonideal_liq, nonideal_gas, Bik
    )

    st.success(f"Converged Bubble Point Temperature: {T_final:.8f} {T_unit_Ant.upper()}")
    st.write("Vapor phase mole fractions yi:")
    st.write(y_final)

