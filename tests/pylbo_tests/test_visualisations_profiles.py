import pylbo


def test_plot_equilibrium_profile(ds_v112):
    pylbo.plot_equilibrium(ds_v112)


def test_plot_continuum_profile(ds_v112):
    pylbo.plot_continua(ds_v112)


def test_plot_equilibrium_balance(ds_v112):
    pylbo.plot_equilibrium_balance(ds_v112)


def test_plot_matrices(ds_v100):
    pylbo.plot_matrices(ds_v100)
