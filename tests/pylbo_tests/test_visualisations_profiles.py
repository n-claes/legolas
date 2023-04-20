import pylbo


def test_plot_equilibrium_profile(ds_v112):
    pylbo.plot_equilibrium(ds_v112)


def test_plot_continuum_profile(ds_v112):
    pylbo.plot_continua(ds_v112)


def test_plot_equilibrium_balance(ds_v112):
    pylbo.plot_equilibrium_balance(ds_v112)


def test_plot_matrices(ds_v100):
    pylbo.plot_matrices(ds_v100)


def test_plot_equilibrium_profile_no_bg(ds_v200_tear_nobg):
    p = pylbo.plot_equilibrium(ds_v200_tear_nobg)
    assert p is None


def test_plot_equilibrium_balance_no_bg(ds_v200_tear_nobg):
    p = pylbo.plot_equilibrium_balance(ds_v200_tear_nobg)
    assert p is None
