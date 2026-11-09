from scipy import signal


def plant(gamma1):
    """Plant transfer function P(s) = 1 / (s + gamma1)."""
    return signal.TransferFunction([1.0], [1.0, gamma1])


if __name__ == "__main__":
    P = plant(gamma1=1.0)
    t, y = signal.step(P)
    print("pole:", P.poles)
    print(f"final value: {y[-1]:.3f} (expect 1/gamma1 = 1)")
