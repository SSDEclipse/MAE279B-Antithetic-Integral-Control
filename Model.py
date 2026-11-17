from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

def closed_loop(gamma1, omega0, Kp, Ki, Ks=1.0, Kf=0.0):
    """Linearized closed-loop transfer function, Eq. (9).
    Kp = 0 -> rAIF, Kf = 0 -> sAIF."""
    num = [omega0 * Kf, omega0 * Ki] if Kf else [omega0 * Ki]
    den = [1.0,
           omega0 + gamma1,
           omega0 * (gamma1 + Kp * Ks),
           omega0 * Ks * Ki]
    return signal.TransferFunction(num, den)


if __name__ == "__main__":
    H = closed_loop(gamma1=1.0, omega0=2.0, Kp=1.0, Ki=2.0)  # sAIF (Kf = 0)
    t, y = signal.step(H, N=1000)
    poles = H.poles
    print("poles:", poles.imag)
    print(f"final value: {y[-1]:.3f} (expect 1/Ks = 1)")

    rl_poles = []
    rl_num_gains = 1000
    rl_kps = np.linspace(0, 1000, rl_num_gains)
    rl_kis = np.linspace(0, 0.5, rl_num_gains)
    for i in range(rl_num_gains):
        H = closed_loop(gamma1=1.0, omega0=1.5, Kp=0, Ki=rl_kis[i])
        rl_poles.append(H.poles)

    plt.figure()
    plt.plot(t, y)
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Output')
    plt.title('Step Response')

    plt.figure()

    for i in range(len(rl_poles[0])):
        reals = []
        imags = []
        for j in range(len(rl_poles)):
            reals.append(rl_poles[j][i].real)
            imags.append(rl_poles[j][i].imag)
        plt.plot(reals, imags, marker='x')
    plt.grid()
    plt.xlabel('Real')
    plt.ylabel('Imaginary')
    plt.title('Root Locus')
    plt.show()



# Note From Luke
# Equation 9 is not the model transfer function. It is the whole system.
# We implement rAIF and sAIF by two different methods:

# rAIF: K_P = 0
# sAIF: K_F = 0

# We should talk about this on thursday. The plant model itslf is very simple.

# Also, I think we missed something about how the paper is set up. They linearize the system
# and reduce it down to one transfer function for analysis purposes

# However, the actual model is a 3 equation system, and we'll need to build this too.
