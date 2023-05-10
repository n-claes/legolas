from __future__ import annotations

import numpy as np
from pylbo.data_containers import LegolasDataSet
from pylbo.utilities.logger import pylboLogger


def get_equilibrium_balance(ds: LegolasDataSet) -> dict:
    balance = {
        "force": _force_balance_1(ds),
        "energy": _energy_balance(ds),
    }
    if _needs_continuity(ds):
        balance["continuity"] = _continuity_balance(ds)
    if _needs_force_balance_2(ds):
        balance["force 2"] = _force_balance_2(ds)
    if _needs_force_balance_3(ds):
        balance["force 3"] = _force_balance_3(ds)
    if _needs_induction_balance_1(ds):
        balance["induction 1"] = _induction_balance_1(ds)
    if _needs_induction_balance_2(ds):
        balance["induction 2"] = _induction_balance_2(ds)
    return balance


def _needs_continuity(ds: LegolasDataSet) -> bool:
    return any(ds.equilibria["v01"] != 0)


def _needs_force_balance_2(ds: LegolasDataSet) -> bool:
    return any(ds.equilibria["v01"] != 0) or any(ds.equilibria["B01"] != 0)


def _needs_force_balance_3(ds: LegolasDataSet) -> bool:
    return _needs_force_balance_2(ds)


def _needs_induction_balance_1(ds: LegolasDataSet) -> bool:
    return _needs_force_balance_2(ds)


def _needs_induction_balance_2(ds: LegolasDataSet) -> bool:
    return _needs_force_balance_2(ds)


def _continuity_balance(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    eps_f = ds.d_scale_factor / ds.scale_factor  # deps / eps
    return (
        bg["drho0"] * bg["v01"]
        + bg["rho0"] * bg["dv01"]
        + bg["rho0"] * bg["v01"] * eps_f
    )


def _force_balance_1(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    eps_f = ds.d_scale_factor / ds.scale_factor
    return (
        bg["drho0"] * bg["T0"]
        + bg["rho0"] * bg["dT0"]
        + bg["B02"] * bg["dB02"]
        + bg["B03"] * bg["dB03"]
        + bg["rho0"] * bg["gravity"]
        - (eps_f) * (bg["rho0"] * bg["v02"] ** 2 - bg["B02"] ** 2)
        + bg["rho0"] * bg["v01"] * bg["dv01"]
    )


def _force_balance_2(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    eps_f = ds.d_scale_factor / ds.scale_factor
    return bg["rho0"] * bg["v01"] * (bg["dv02"] + bg["v02"] * eps_f) - bg["B01"] * (
        bg["dB02"] + bg["B02"] * eps_f
    )


def _force_balance_3(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    return bg["rho0"] * bg["v01"] * bg["dv03"] - bg["B01"] * bg["dB03"]


def _energy_balance(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    eps = ds.scale_factor
    deps = ds.d_scale_factor
    dkappa_perp_dr = bg.get("dkappa_perp_dr", None)
    if dkappa_perp_dr is None:
        dkappa_perp_dr = _derivative_from_gradient(
            fname="dkappa_perp_dr",
            fname_prim="kappa_perp",
            bg=bg,
            with_respect_to=ds.grid_gauss,
        )
    Kp = _get_conduction_prefactor(ds)
    dKp = _get_conduction_prefactor_derivative(ds)
    return (
        bg["T0"] * bg["rho0"] * (deps * bg["v01"] + eps * bg["dv01"]) / eps
        + bg["rho0"] * bg["L0"]
        - bg["B01"] ** 2 * (Kp * bg["dT0"] + dKp * bg["T0"])
        - (1 / eps)
        * (
            deps * bg["kappa_perp"] * bg["dT0"]
            + eps * dkappa_perp_dr * bg["dT0"]
            + eps * bg["kappa_perp"] * bg["ddT0"]
        )
        + (1 / (ds.gamma - 1)) * bg["dT0"] * bg["rho0"] * bg["v01"]
    )


def _induction_balance_1(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    return bg["B01"] * bg["dv02"] - bg["B02"] * bg["dv01"] - bg["dB02"] * bg["v01"]


def _induction_balance_2(ds: LegolasDataSet) -> np.ndarray:
    bg = ds.equilibria
    eps_f = ds.d_scale_factor / ds.scale_factor
    return (
        bg["B01"] * bg["dv03"]
        - bg["dB03"] * bg["v01"]
        - bg["B03"] * bg["dv01"]
        + eps_f * (bg["B01"] * bg["v03"] - bg["B03"] * bg["v01"])
    )


def _get_conduction_prefactor(ds: LegolasDataSet) -> np.ndarray:
    if not _ds_is_mhd(ds):
        return np.zeros_like(ds.grid_gauss)
    bg = ds.equilibria
    return (bg["kappa_para"] - bg["kappa_perp"]) / bg["B0"] ** 2


def _get_conduction_prefactor_derivative(ds: LegolasDataSet) -> np.ndarray:
    if not _ds_is_mhd(ds):
        return np.zeros_like(ds.grid_gauss)
    bg = ds.equilibria
    dB0 = (bg["B02"] * bg["dB02"] + bg["B03"] * bg["dB03"]) / bg["B0"]
    dkappa_para_dr = bg.get("dkappa_para_dr", None)
    if dkappa_para_dr is None:
        dkappa_para_dr = _derivative_from_gradient(
            fname="dkappa_para_dr",
            fname_prim="kappa_para",
            bg=bg,
            with_respect_to=ds.grid_gauss,
        )
    dkappa_perp_dr = bg.get("dkappa_perp_dr", None)
    if dkappa_perp_dr is None:
        dkappa_perp_dr = _derivative_from_gradient(
            fname="dkappa_perp_dr",
            fname_prim="kappa_perp",
            bg=bg,
            with_respect_to=ds.grid_gauss,
        )
    return (
        (dkappa_para_dr - dkappa_perp_dr) * bg["B0"]
        - 2 * (bg["kappa_para"] - bg["kappa_perp"]) * dB0
    ) / bg["B0"] ** 3


def _derivative_from_gradient(
    fname: str, fname_prim: str, bg: dict, with_respect_to: np.ndarray
) -> np.ndarray:
    pylboLogger.warning(
        f"field with name '{fname}' not found in datfile. "
        f"Deriving numerically from '{fname_prim}' using np.gradient."
    )
    return np.gradient(bg[fname_prim], with_respect_to, edge_order=2)


def _ds_is_mhd(ds: LegolasDataSet) -> bool:
    return "mhd" in ds.header.get("physics_type", None)
