---
title: Equations & physics
layout: single
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "calculator"
last_modified_at: 2020-10-28
---

In general, Legolas solves the linearised MHD equations using a Fourier analysis in the ignorable
coordinates. The full system of equations is given below, including flow, gravity, resistivity,
optically thin radiative losses and anisotropic thermal conduction.

For a detailed discussion we refer to our [code paper](https://arxiv.org/abs/2010.14148).

## System of equations

$$
\newcommand{\bfv}{\boldsymbol{v}}
\newcommand{\bfb}{\boldsymbol{B}}
\newcommand{\bfg}{\boldsymbol{g}}
\newcommand{\bfkappa}{\boldsymbol{\kappa}}
\newcommand{\HL}{\mathscr{L}}
\begin{align}
	\frac{\partial \rho}{\partial t} &= -\nabla \cdot (\rho\bfv) \\
	\rho\frac{\partial \bfv}{\partial t} &= -\nabla p - \rho\bfv \cdot \nabla \bfv + (\nabla \times \bfb) \times \bfb + \rho\bfg	\\
	\rho\frac{\partial T}{\partial t} &= -\rho\bfv \cdot \nabla T - (\gamma - 1)p\nabla\cdot\bfv - (\gamma - 1)\rho\HL + (\gamma - 1)\nabla\cdot(\bfkappa \cdot \nabla T) + (\gamma - 1)\eta(\nabla \times \bfb)^2\\
	\frac{\partial \bfb}{\partial t} &= \nabla \times (\bfv \times \bfb) - \nabla \times (\eta\nabla \times \bfb)
\end{align}
$$

## Equilibrium state
Legolas requires a time-independent background equilibrium of the form

$$
\newcommand{\bey}{\boldsymbol{u}_2}
\newcommand{\bez}{\boldsymbol{u}_3}
\begin{aligned}
    \rho_0 &= \rho_0(u_1)		\\
    p_0 &= p_0(u_1) 			\\
    T_0 &= T_0(u_1) 			\\
\end{aligned}
\qquad
\begin{aligned}
    \bfv_0 &= v_{02}(u_1)\bey + v_{03}(u_1)\bez	\\
    \bfb_0 &= B_{02}(u_1)\bey + B_{03}(u_1)\bez \\
\end{aligned}
$$

where $(u_1, u_2, u_3)$ denotes the coordinate system: $(x, y, z)$ in Cartesian and $(r, \theta, z)$ in cylindrical
geometries.

## Equilibrium requirements
By extension, the combination of the equilibrium state and the full system of equations yields conditions that have
to be fulfilled. In particular, every equilibrium state must satisfy the following two equations

$$
\newcommand{\eps}{\varepsilon}
\begin{align}
	\left(p_0 + \frac{1}{2}\bfb_0^2\right)' + \rho_0 g - \frac{\eps'}{\eps}\left(\rho_0v_{02}^2 - B_{02}^2\right) &= 0 \\
	\frac{1}{\eps}\left(\eps \kappa_\bot T_0'\right)' - \rho_0\HL_0 &= 0
\end{align}
$$

where $\varepsilon = r$ and $\varepsilon' = 1$ in cylindrical geometries, and $\varepsilon = 1$ and $\varepsilon' = 0$
in Cartesian slabs. The prime denotes the derivative with respect to $u_1$.

Legolas explicitly checks these conditions at runtime and will raise a warning if these are not fulfilled,
giving you additional information on how large the discrepancies are and where they occur with respect to the grid.

## Units
All equations in Legolas are in dimensionless form, as is common practice when dealing with (M)HD.
We adopt a standard reference value of 2 for plasma beta, such that we can write $\beta/2 = \mu p / B^2 = 1$ with
$\mu$ the magnetic constant, equal to $4\pi$ in cgs units.
A unit magnetic field $B_{unit}$ and unit length $L_{unit}$ should _always_ be specified, yielding a unit pressure of

$$
p_{unit} = \frac{B_{unit}^2}{\mu}
$$

Then, specifying a unit density $\rho_{unit}$ fixes the unit temperature $T_{unit}$ (or vice-versa)
through the ideal gas law

$$
p_{unit} = \mathcal{R}_{specific}T_{unit}\rho_{unit}
$$

Hence we only have **three** degrees of freedom: the unit magnetic field $B_{unit}$ and unit length $L_{unit}$
should be set, together with _either_ the unit density $\rho_{unit}$ OR unit temperature $T_{unit}$.
All other unit normalisations follow from the three degrees of freedom and are given by

$$
\begin{align}
    n_{unit} &= \frac{\rho_{unit}}{m_p} \qquad\qquad v_{unit} = \frac{B_{unit}}{\sqrt{\mu\rho_{unit}}}    \\
    t_{unit} &= \frac{L_{unit}}{v_{unit}} \qquad\qquad    \mathscr{L}_{unit} = \frac{p_{unit}}{t_{unit}n_{unit}^2}    \\
    \kappa_{unit} &= \frac{\rho_{unit}L_{unit}v_{unit}^3}{T_{unit}}
\end{align}
$$

which, from left to right, read as unit numberdensity, unit velocity, unit time, unit luminosity and unit
conduction.

**Note:** The unit normalisations are only used when radiative cooling, thermal conduction or temperature-dependent
resistivity is included. We always set base values though (as one should), which are given in cgs units by
$B_{unit} = 10$ Gauss, $L_{unit} = 10^9$ cm and $T_{unit} = 10^6$ K.
{: .notice--info}
