import numpy as np

def exact_solution(x, t):
    """Calculates the exact solution of the advection equation."""
    return np.exp(-20.0 * ((x - t) - 2.0)**2) + np.exp(-((x - t) - 5.0)**2)

def solve_ftcs(u_initial, m, delta_x, delta_t, t_final):
    """Forward-Time, Central-Space (Unstable for pure advection)."""
    U_fe = u_initial.copy()
    U_new = np.empty(U_fe.shape)
    t = 0
    while t < t_final:
        dt = min(delta_t, t_final - t)
        nu = dt / delta_x
        for i in range(m):
            left = U_fe[-1] if i == 0 else U_fe[i-1]
            right = U_fe[0] if i == m-1 else U_fe[i+1]
            U_new[i] = U_fe[i] - 0.5 * nu * (right - left)
        U_fe[:] = U_new
        t += dt
    return U_fe

def solve_lax_friedrichs(u_initial, m, delta_x, delta_t, t_final):
    """Lax-Friedrichs Scheme."""
    U_lf = u_initial.copy()
    U_new = np.empty(U_lf.shape)
    t = 0
    while t < t_final:
        dt = min(delta_t, t_final - t)
        nu = dt / delta_x
        for i in range(m):
            left = U_lf[-1] if i == 0 else U_lf[i-1]
            right = U_lf[0] if i == m-1 else U_lf[i+1]
            U_new[i] = 0.5 * (left + right) - 0.5 * nu * (right - left)
        U_lf[:] = U_new
        t += dt
    return U_lf

def solve_upwind(u_initial, m, delta_x, delta_t, t_final):
    """First-Order Upwind Scheme."""
    U_up = u_initial.copy()
    U_new = np.empty(U_up.shape)
    t = 0
    while t < t_final:
        dt = min(delta_t, t_final - t)
        nu = dt / delta_x
        for i in range(m):
            left = U_up[-1] if i == 0 else U_up[i-1]
            U_new[i] = U_up[i] - nu * (U_up[i] - left)
        U_up[:] = U_new
        t += dt
    return U_up

def solve_lax_wendroff(u_initial, m, delta_x, delta_t, t_final):
    """Lax-Wendroff Scheme (2nd Order)."""
    U_lw = u_initial.copy()
    U_new = np.empty(U_lw.shape)
    t = 0
    while t < t_final:
        dt = min(delta_t, t_final - t)
        nu = dt / delta_x
        for i in range(m):
            left  = U_lw[i-1] if i > 0 else U_lw[-1]
            right = U_lw[i+1] if i < m-1 else U_lw[0]
            U_new[i] = (U_lw[i] - 0.5 * nu * (right - left) + 
                        0.5 * nu**2 * (right - 2*U_lw[i] + left))
        U_lw[:] = U_new
        t += dt
    return U_lw

def solve_beam_warming(u_initial, m, delta_x, delta_t, t_final):
    """Beam-Warming Scheme (2nd Order Upwind)."""
    U_bw = u_initial.copy()
    U_new = np.empty(U_bw.shape)
    t = 0
    while t < t_final:
        dt = min(delta_t, t_final - t)
        nu = dt / delta_x
        for i in range(m):
            U_i   = U_bw[i]
            U_im1 = U_bw[i-1] if i > 0 else U_bw[-1]
            U_im2 = U_bw[i-2] if i > 1 else U_bw[-2 + i]
            U_new[i] = (U_i - 0.5 * nu * (3*U_i - 4*U_im1 + U_im2) + 
                        0.5 * nu**2 * (U_i - 2*U_im1 + U_im2))
        U_bw[:] = U_new
        t += dt
    return U_bw
