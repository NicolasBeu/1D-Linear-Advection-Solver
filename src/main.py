import numpy as np
import matplotlib.pyplot as plt
import os
from src.solvers import (exact_solution, solve_ftcs, solve_lax_friedrichs, 
                         solve_upwind, solve_lax_wendroff, solve_beam_warming)

def main():
    # Setup numerical parameters
    m = 501
    x = np.linspace(0.0, 25.0, m)
    delta_x = 25.0 / (m - 1)
    cfl = 0.8
    delta_t = cfl * delta_x
    t_final = 17.0

    # Initial Condition
    U_initial = exact_solution(x, 0)

    # Run Solvers
    print("Running numerical simulations...")
    U_fe = solve_ftcs(U_initial, m, delta_x, delta_t, t_final)
    U_lf = solve_lax_friedrichs(U_initial, m, delta_x, delta_t, t_final)
    U_up = solve_upwind(U_initial, m, delta_x, delta_t, t_final)
    U_lw = solve_lax_wendroff(U_initial, m, delta_x, delta_t, t_final)
    U_bw = solve_beam_warming(U_initial, m, delta_x, delta_t, t_final)

    # Visualization
    fig, axes = plt.subplots(figsize=(12,6))
    axes.plot(x, exact_solution(x, t_final), 'k', linewidth=2, label='Exact Solution')
    axes.plot(x, U_fe, 'r.', markersize=3, label='FTCS')
    axes.plot(x, U_lf, 'b--', linewidth=1, label='Lax-Friedrichs')
    axes.plot(x, U_up, 'g.', markersize=3, label='Upwind')
    axes.plot(x, U_lw, 'm-.', linewidth=1, label='Lax-Wendroff')
    axes.plot(x, U_bw, 'c:', linewidth=1, label='Beam-Warming')

    axes.set_xlim((15.0, 25.0))
    axes.set_ylim((-0.1, 1.1))
    axes.set_xlabel('x')
    axes.set_ylabel('u(x,t)')
    axes.set_title(f'Numerical Methods Comparison at t = {t_final}')
    axes.legend()

    # Save plot to assets folder
    os.makedirs('assets', exist_ok=True)
    plt.tight_layout()
    plt.savefig('assets/methods_comparison.png', dpi=300)
    print("Simulation complete. Plot saved to assets/methods_comparison.png")

if __name__ == "__main__":
    main()
