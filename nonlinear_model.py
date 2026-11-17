# Coupled ODE solver
# Luke
# Most code written by hand. Claude did a lot of debugging, helped finding a lot of the paramaters,
# and implemented some of the simulation setup stuff. 

# Paramaters may not all be provided in the paper.
# I believe there is also a github repo with all the code for this project in matlab, so we could definitely 
# scrape some of the info from there

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def nonlinear_model(t, state, model_type,  eta, theta, mu, gamma1, alpha=None, k=None, n=None):

    x, z1, z2 = state


    if model_type == "sAIF":
        u = alpha / (1+ (z2/k)**n)
    elif model_type == "rAIF":
        u = k * z1
    else:
        raise ValueError(f"model_type must be rAIF or sAIF")
    
   
    x_dot = u - gamma1 * x
    z1_dot = mu - eta * z1 * z2
    z2_dot = theta * x - eta * z1 * z2

    

    return [x_dot, z1_dot, z2_dot]



# It's also not totally clear if they use the linearized model or the full ODE setup. 


model_type = "rAIF"
gamma1 = 1.0
theta = 1.0
mu0 = 5.0
mu1 = 2 * mu0
kappa_x = 1e-5
n = 2 # This is sometimes varied!

u_bar = gamma1 * mu1 / theta
u_bar0 = gamma1 * mu0 / theta

# Some values are computed:
if model_type == "rAIF":
    Wo = 1.5
    KI = 0.175
elif model_type == "sAIF":
    a = 1.0
    Ks = theta
    Wo = 3 * a - gamma1
    KI = a**3 / (Ks * Wo)
    KP = (3 * a**2 - gamma1 * Wo) / (Ks * Wo)
    if not (KP < n * u_bar / mu1 and KI < Wo * KP * (1 - mu1 * KP / (n * u_bar))):
        raise ValueError("Gains outside repression coverage (Eq. 8): try a larger n or a different a")


if model_type == "rAIF":
    # From Equation 24
    eta = Wo * KI / u_bar
    k = u_bar / mu1 * Wo / 2 * (1 - (1 - 4 * mu1 / u_bar * KI / Wo)**0.5)
    # Interestingly, the paper doesn't say which to use. Could be something to comment on in the presentation.
    alpha = None

    z1_0 = u_bar0 / k
    z2_0 = mu0 / (eta * z1_0)

elif model_type == "sAIF":
    # From Equation 30
    r = KI / KP
    eta = 1 / mu1 * r * (Wo - r)
    alpha = u_bar / (1 - mu1 / (n * u_bar) * Wo * KP / (Wo - r))
    kappa = mu1 / (Wo - r) * (n * u_bar / mu1 * (Wo - r) / (Wo * KP) - 1)**(1 / n)
    k = kappa

    z2_0 = kappa * (alpha / u_bar0 - 1)**(1 / n)
    z1_0 = mu0 / (eta * z2_0)

x0 = mu0 / theta
print(f"{model_type}: Wo={Wo}, KI={KI}, eta={eta:.4g}, k={k:.4g}, alpha={alpha}")

sol = solve_ivp(nonlinear_model, (0, 20), [x0, z1_0, z2_0],
                args=(model_type, eta, theta, mu1, gamma1, alpha, k, n),
                method="LSODA", t_eval=np.linspace(0, 20, 500))

print(f"final x = {sol.y[0, -1]:.3f} (expect mu1/theta = {mu1 / theta})")
plt.figure()
plt.plot([-3, 0], [x0, x0], color="C0")
plt.plot(sol.t, sol.y[0], color="C0", label="x (output)")
plt.axhline(mu1 / theta, color="k", linestyle="--", label="setpoint")
plt.xlabel("Time")
plt.ylabel("Output")
plt.legend()
plt.grid()

plt.figure()
plt.plot(sol.t, sol.y[1], label="z1")
plt.plot(sol.t, sol.y[2], label="z2")
plt.xlabel("Time")
plt.legend()
plt.grid()
plt.show()
