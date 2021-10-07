---
title: Equations
layout: single
classes: wide
sidebar:
  nav: "leftcontents"
toc: true
toc_icon: "chevron-circle-down"
last_modified_at: 2021-07-27
---

On this page we give a small overview of the system of equations solved by Legolas. We use the
full set of MHD equations, linearised around a dynamic background (that is, including flow).
Physical effects include flow, (external) gravity, resistivity, optically thin radiative losses,
anisotropic thermal conduction, viscosity and Hall effects.

## System of equations
<!-- mathematical symbols used on this page (spaces needed before and after $$ to render properly) -->

$$
\newcommand{\vbf}{\mathbf{v}}
\newcommand{\gbf}{\mathbf{g}}
\newcommand{\bbf}{\mathbf{B}}
\newcommand{\HL}{\mathscr{L}}
\newcommand{\HH}{\mathcal{H}}
\newcommand{\kappabf}{\boldsymbol{\kappa}}
\newcommand{\unit}[1]{\mathbf{e}_{#1}}

\begin{gather}
	\frac{\partial \rho}{\partial t} = -\nabla \cdot (\rho \vbf), \\
	\rho\frac{\partial \vbf}{\partial t} = -\nabla p - \rho \vbf \cdot \nabla \vbf + (\nabla \times \bbf) \times \bbf
                    + \rho\gbf + \mu\left[\nabla^2\vbf + \frac{1}{3}\nabla(\nabla \cdot \vbf)\right], \\
	\begin{aligned}
		\rho\frac{\partial T}{\partial t} = &-\rho \vbf\cdot\nabla T - (\gamma - 1)p\nabla \cdot \vbf- (\gamma - 1)\rho\HL \\
                &+ (\gamma - 1)\nabla \cdot (\kappabf \cdot \nabla T)
				+ (\gamma - 1)\eta(\nabla \times \bbf)^2 + \mu\left|\nabla \vbf \right|^2,
	\end{aligned} \\
	\begin{aligned}
		\frac{\partial \bbf}{\partial t} = ~&\nabla \times (\vbf \times \bbf) - \nabla \times (\eta\nabla \times \bbf) \\
				&-\nabla\times \left[ \frac{\eta_H}{\rho_e}(\nabla \times \bbf) \times \bbf - \frac{\eta_H}{\rho_e}\nabla p
                                        + \frac{\eta_e}{\rho_e}\frac{\partial (\nabla \times \bbf)}{\partial t} \right].
	\end{aligned}
\end{gather}
$$

## Physical effects
Below we briefly describe the various physical effects that can be added. For more information on how to add these in your particular
setup take a look at the [user submodule](../../general/own_setup#including-additional-physics) page.

### External gravity
An external gravitational field can be added, which can either be constant or dependent on position. In Cartesian geometry we assume this to be
$x$-dependent and aligned with height:

$$
\gbf = -g(x)\unit{x}.
$$

A similar profile can be used in cylindrical geometries, e.g. for accretion disk setups.

### Resistivity
Resistivity can either be constant over the entire grid, or a more general temperature-dependent profile given
by the Spitzer resistivity:

$$
\eta = \frac{4\sqrt{2\pi}}{3}\frac{Z_\text{ion}e^2\sqrt{m_e}\ln(\lambda)}{(4\pi\epsilon_0)^2(k_B T)^{3/2}},
$$

with the Coulomb logarithm $\ln(\lambda) \approx 22$.

### Optically thin radiative losses
Radiative cooling is governed by the heat-loss function, specified as the difference between energy gains and energy losses

$$
\HL = \rho\Lambda(T) - \HH,
$$

The function $\Lambda(T)$ here is called the _cooling curve_, which is a tabulated set of values resulting from detailed molecular calculations.
Legolas has multiple cooling curves from existing literature implemented which are interpolated at high resolution and sampled on the given grid.
We have also included analytical, piecewise prescriptions (e.g. the `rosner` curve). See the [physicslist](../../general/parameter_file/#physicslist) for an overview of the different options.

The function $\HH$ is called the _heating term_, and it is currently not (yet) possible to specify this. We assume that this term only depends on the equilibrium background, meaning that it is
constant in time but possibly varying in space and as such that it balances out the cooling contribution to ensure thermal equilibrium.

### Thermal conduction
The thermal conduction prescription relies on a tensor representation to model the anisotropy in MHD and is given by

$$
\kappabf = \kappa_\parallel\unit{B}\unit{B} + \kappa_\bot(\boldsymbol{I} - \unit{B}\unit{B}),
$$

with $\boldsymbol{I}$ the unit tensor and $\unit{B} = \bbf / B$ a unit vector along the magnetic field. For the actual coefficients the Spitzer conductivity is assumed:

$$
\begin{align}
\kappa_\parallel &\approx 8 \times 10^{-7}T^{5/2}~\text{erg cm$^{-1}$s$^{-1}$K$^{-1}$},	\\
\kappa_\bot &\approx 4 \times 10^{-10} n^2B^{-2}T^{-3}\kappa_\parallel,
\end{align}
$$

Alternatively, both of these can be independently set to constant values.

### Viscosity
The viscosity addition to the momentum equation is usually given in a full tensor form, but in Legolas this is simplified to good approximation to

$$
\mathbf{F}_{\mathrm{v}} = -\nabla\cdot\boldsymbol{\Pi} \simeq \mu \left[ \nabla^2\vbf + \frac{1}{3} \nabla\left( \nabla\cdot\vbf \right) \right],
$$

where the constant $\mu$ denotes the dynamic viscosity. The viscous heating term in the energy equation is expanded using a Frobenius norm and written as

$$
\left|\nabla \vbf\right|^2 = \sum_{i=1}^{3}\sum_{j=1}^{3}\left(\nabla \vbf\right)_{ij}^2
$$

The heating term can be toggled on or off independently of the main viscosity addition in the momentum equation,
see the [physicslist](../../general/parameter_file/#physicslist) for details on how to do this.


### Hall MHD
Legolas allows for the addition of Hall effects, including electron contributions:

$$
\underbrace{\frac{\eta_H}{\rho}(\nabla \times \bbf) \times \bbf}_\text{main Hall term}, \qquad\qquad
\underbrace{\frac{\eta_H}{\rho}\nabla p}_\text{electron pressure}, \qquad\qquad
\underbrace{\frac{\eta_e}{\rho}\frac{\partial (\nabla \times \bbf)}{\partial t}}_\text{electron inertia}.
$$

The main Hall contribution, electron pressure term and electron inertia term are treated as separate entities, meaning that
these terms can be toggled on or off independently of one another. See the [physicslist](../../general/parameter_file/#physicslist) for details.
