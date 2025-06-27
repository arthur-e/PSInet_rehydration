Sundown until refilling (predawn) on the X-axis
Average water potential in that period (or Water potential at rehydration?) on the Y axis

Keff = -log(psi_t - psi_pd / (psi_gs - psi_pd)) (C / t)

- C is the average leaf capacitance
- psi_gs is the psi_leaf at nighttime stomatal closure (sundown)
- psi_pd at predawn or highest leaf water potential

Rehydration time constant is Keff/C

tau = R_0 * C
tau = C / Keff

In the case of an electrical circuit (RC circuit)...

$$
V(t) = V_0 (1 - e^{-t/\tau})\\
\tau = RC
$$

Translating this to plant hydraulics:

$$
\Psi(t) = \Psi_{soil}(1 - e^{-t/\tau})\\
\tau = \frac{C}{K}
$$

We hope that $\Psi_{soil}$ can be represented by the pre-dawn plant water potential, in the absence of nighttime transpiration. If there is nocturnal transpiration:
$$
\Psi_{pd} = \Psi_{soil} - \frac{E_n}{K_n}
$$

Where $E_n$ is the nighttime transpiration and $K_n$ is the nighttime hydraulic conductance.

Giving us:
$$
\Psi(t) = (\Psi_{soil} - \frac{E_n}{K_n})(1 - e^{(t\times K_n)/ C})
$$

Where $E_n$ might be replaced with the product of VPD and stomatal conductance.

NOTE: Should be one plus the exponential when $\Psi_{soil}$ is negative (as water potential is).


References
==========

Donovan et al. (2001, Oecologia), "Predawn plant water potential does not necessarily..."
