---
title: Unit normalisations
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
last_modified_at: 2023-04-23
---

All equations in Legolas are in dimensionless form, as is common practice when dealing with (M)HD.
As usual we have **three** degrees of freedom.

## Mean molecular weight
Unit normalisations depend on the molecular weight $\bar{\mu}$, and in Legolas we usually distinguish between two cases:

- **Electron-proton plasma**: $\bar{\mu} = 0.5$, this is the default case.

    $$
    \bar{\mu} = \dfrac{m_e n_e + m_i n_i}{n_e + n_i} \simeq \dfrac{m_i n_i}{n_e + n_i}
              = \dfrac{1}{2}m_i \rightarrow \bar{\mu} = \dfrac{1}{2}
    $$

- **Pure proton plasma**: $\bar{\mu} = 1$, in this case the molecular weight should be explicitly set.

    $$
    \bar{\mu} = \dfrac{m_i n_i}{n_i} = m_i \rightarrow \bar{\mu} = 1
    $$

## Normalisations
Legolas has three options to specify units, all in cgs. In what follows $m_p$ denotes the proton mass,
$k_B$ the Boltzmann constant, and $\mu_0 = 4\pi$ the magnetic constant.

1. Reference unit density, unit magnetic field and unit length $(\rho_u, B_u, L_u)$, then

   $$
   p_u = \frac{B_u^2}{\mu_0}, \quad
   T_u = \frac{\bar{\mu} p_u m_p}{k_B \rho_u}, \quad
   n_u = \frac{\rho_u}{m_p}, \quad
   v_u = \frac{B_u}{\sqrt{\mu_0 \rho_u}}.
   $$

2. Reference unit temperature, unit magnetic field and unit length $(T_u, B_u, L_u)$, then

   $$
   p_u = \frac{B_u^2}{\mu_0}, \quad
   \rho_u = \frac{\bar{\mu} p_u m_p}{k_B T_u}, \quad
   n_u = \frac{\rho_u}{m_p}, \quad
   v_u = \frac{B_u}{\sqrt{\mu_0 \rho_u}}.
   $$

3. Reference unit numberdensity, unit temperature and unit length $(n_u, T_u, L_u)$, then

   $$
   p_u = \bar{\mu} n_u k_B T_u, \quad
   \rho_u = m_p n_u, \quad
   v_u = \sqrt{\frac{p_u}{\rho_u}}, \quad
   B_u = \sqrt{\mu_0 p_u}.
   $$

All other normalisations follow from those above and are given by
- unit mass: $M_u = \rho_u L_u^3$
- unit time: $t_u = \dfrac{L_u}{v_u}$
- unit resistivity: $\eta_u = \dfrac{L_u^2}{t_u}$
- unit cooling curve: $\Lambda_u = \dfrac{p_u}{t_u n_u^2}$
- unit conduction: $\kappa_u = \dfrac{\rho_u L_u v_u^3}{T_u}$


<i class="fas fa-lightbulb" aria-hidden="true"></i>
**Note:** the unit normalisations are only relevant when radiative cooling, thermal conduction or temperature-dependent resistivity is included.
We always set base values though (as one should), which are set using option 2. with default values
$B_u = 10$ G, $L_u = 10^9$ cm and $T_u = 10^6$ K.
{: .notice--info}
