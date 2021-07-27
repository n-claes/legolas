---
title: Unit normalisations
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2021-07-27
---

All equations in Legolas are in dimensionless form, as is common practice when dealing with (M)HD.

As usual we have **three** degrees of freedom: the unit magnetic field $B_u$ and unit length $L_u$ should be set,
together with _either_ the unit density $\rho_u$ **OR** unit temperature $T_u$.

## Mean molecular weight
Unit normalisations depend on the molecular weight $\bar{\mu}$, and in Legolas we usually distinguish between two cases:

- **Pure proton plasma**: $\bar{\mu} = 1$, this is the default case.

    $$
    \bar{\mu} = \dfrac{m_i n_i}{n_i} = m_i \rightarrow \bar{\mu} = 1
    $$

- **Electron-proton plasma**: $\bar{\mu} = 0.5$, in this case the molecular weight should be explicitly set.

    $$
    \bar{\mu} = \dfrac{m_e n_e + m_i n_i}{n_e + n_i} \simeq \dfrac{m_i n_i}{n_e + n_i}
              = \dfrac{1}{2}m_i \rightarrow \bar{\mu} = \dfrac{1}{2}
    $$

## Normalisations
We have two options: either the unit density is set or the unit temperature is set. In both cases a unit magneticfield
and a unit length should be set as well, so the unit pressure is given by

$$
p_u = \dfrac{B_u^2}{\mu_0},
$$

taking a plasma-beta reference value of 2. The magnetic constant is given by $\mu_0$ ($= 4\pi$ in cgs).

We then have the following options:

- If $\rho_u$ chosen along with $B_u$ and $L_u$:

    $$
    T_u = \dfrac{\bar{\mu}m_p p_u}{k_B \rho_u},
    $$

- If $T_u$ chosen along with $B_u$ and $L_u$:

    $$
    \rho_u = \dfrac{\bar{\mu}m_p p_u}{k_B T_u},
    $$

All other normalisations follow from these, where the Alfvén speed is assumed as reference velocity:

$$
\begin{gather}
    m_u = \rho_u L_u, \\
    n_u = \dfrac{\rho_u}{m_p}, \\
    v_u = \dfrac{B_u}{\sqrt{\mu_0 \rho_u}}, \\
    t_u = \dfrac{L_u}{v_u}, \\
    \Lambda_u = \dfrac{p}{t_u n_u^2}, \\
    \kappa_u = \dfrac{\rho_u L_u v_u^3}{T_u},
\end{gather}
$$

which, from top to bottom, read as mass unit, numberdensity unit, velocity unit, time unit, cooling curve unit and
conduction unit. The proton mass is given by $m_p$ and $k_B$ denotes the Boltzmann constant.

<i class="fas fa-lightbulb" aria-hidden="true"></i>
**Note:** The unit normalisations are only relevant when radiative cooling, thermal conduction or temperature-dependent
resistivity is included. We always set base values though (as one should), which are given in cgs units by
$B_u = 10$ G, $L_u = 10^9$ cm and $T_u = 10^6$ K. If no normalisations are specified these are used by default.
{: .notice--info}