import streamlit as st
import math

st.title("Bubble Point Calculation App")

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

n = st.sidebar.number_input("Number of components", min_value=1, value=2)

nonideal_liq = st.sidebar.radio("Liquid phase non-ideal?", ["Yes", "No"]) == "Yes"
nonideal_gas = st.sidebar.radio("Gas phase non-ideal?", ["Yes", "No"]) == "Yes"

# Wilson parameters
a = [[0.0]*n for _ in range(n)]
V = [1.0]*n

if nonideal_liq:
    st.sidebar.subheader("Wilson Parameters")
    for i in range(n):
        for j in range(n):
            if i != j:
                a[i][j] = st.sidebar.number_input(f"a{i+1}{j+1}", value=0.0, format="%.8f")
    V = [st.sidebar.number_input(f"V{i+1}", value=1.0, format="%.8f") for i in range(n)]

# Virial coefficients
Bik = [[0.0]*n for _ in range(n)]
if nonideal_gas:
    st.sidebar.subheader("Virial Coefficients")
    for i in range(n):
        for j in range(i, n):
            val = st.sidebar.number_input(f"B{i+1}{j+1}", value=0.0, format="%.8f")
            Bik[i][j] = Bik[j][i] = val

A = [st.sidebar.number_input(f"A{i+1}", value=8.07131, format="%.8f") for i in range(n)]
B = [st.sidebar.number_input(f"B{i+1}", value=1730.63, format="%.8f") for i in range(n)]
C = [st.sidebar.number_input(f"C{i+1}", value=233.426, format="%.8f") for i in range(n)]

Ptot = st.sidebar.number_input("System Pressure", value=101.325, format="%.8f")
R = st.sidebar.number_input("Gas Constant R", value=8.314, format="%.8f")

T_unit = st.sidebar.selectbox("Temperature unit", ["C", "K"]).lower()
P_unit = st.sidebar.selectbox("Pressure unit", ["mmHg", "bar", "kPa"]).lower()

T_unit_Ant = T_unit
P_unit_Ant = P_unit


P = convert_pressure_to_kpa(Ptot, P_unit)

k_index = st.sidebar.number_input("Fixed component", min_value=1, max_value=n, value=1) - 1

x = [st.sidebar.number_input(f"x{i+1}", value=1/n, format="%.8f") for i in range(n)]

tol = st.sidebar.number_input("Temperature tolerance (%)", value=0.01, format="%.8f")

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

