import numpy as np
import pylbo
import pytest

nb_u2_pts = 10
nb_u3_pts = 10

u2cart = np.linspace(0, 1, nb_u2_pts)
u2cyl = np.linspace(0, 2 * np.pi, nb_u2_pts)
u3 = np.linspace(0, 3, nb_u3_pts)
omega_cart = 1.19029136 + 3.75969744j
omega_cyl = 0.01746995 + 0.02195201j


@pytest.fixture(scope="function")
def vtkdata_cart(ds_v121_rti_khi):
    return pylbo.prepare_vtk_export(ds_v121_rti_khi, omega_cart, u2=u2cart, u3=u3)


@pytest.fixture(scope="function")
def vtkdata_cyl(ds_v121_magth):
    return pylbo.prepare_vtk_export(ds_v121_magth, omega_cart, u2=u2cyl, u3=u3)


@pytest.fixture(scope="function")
def u1cart(ds_v121_rti_khi):
    return ds_v121_rti_khi.ef_grid


@pytest.fixture(scope="function")
def u1cyl(ds_v121_magth):
    return ds_v121_magth.ef_grid


def test_vtk_data_cart(vtkdata_cart, u1cart):
    assert vtkdata_cart.data is not None
    assert vtkdata_cart.dims == (len(u1cart), nb_u2_pts, nb_u3_pts)


def test_vtk_data_cyl(vtkdata_cyl, u1cyl):
    assert vtkdata_cyl.data is not None
    assert vtkdata_cyl.dims == (len(u1cyl), nb_u2_pts, nb_u3_pts)


def test_vtk_data_cart_coords(vtkdata_cart, u1cart):
    u1, u2, z = vtkdata_cart.get_coordinate_data()
    dims = (len(u1cart), nb_u2_pts, nb_u3_pts)

    u1_expect, u2_expect, z_expect = np.meshgrid(u1cart, u2cart, u3, indexing="ij")
    for u in (u1, u2, z):
        assert u.shape == dims
    assert np.allclose(u1, u1_expect)
    assert np.allclose(u2, u2_expect)
    assert np.allclose(z, z_expect)


def test_vtk_data_cyl_coords(vtkdata_cyl, u1cyl):
    u1, u2, z = vtkdata_cyl.get_coordinate_data()
    dims = (len(u1cyl), nb_u2_pts, nb_u3_pts)

    u1_expect, u2_expect, z_expect = np.meshgrid(u1cyl, u2cyl, u3, indexing="ij")
    for u in (u1, u2, z):
        assert u.shape == dims
    assert np.allclose(u1, u1_expect * np.cos(u2_expect))
    assert np.allclose(u2, u1_expect * np.sin(u2_expect))
    assert np.allclose(z, z_expect)


def test_vtk_write_cart_rho(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_rho.vtk"
    vtkdata_cart.export_to_vtk(filename=file, time=0, names="rho")
    expected = tmpdir / "test_cart_rho_t0000.vtk"
    assert expected.is_file()


def test_vtk_write_cart_rho_bg(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_rho_bg.vtk"
    vtkdata_cart.export_to_vtk(filename=file, time=0, names="rho", bg_names=["rho0"])
    expected = tmpdir / "test_cart_rho_bg_t0000.vtk"
    assert expected.is_file()


def test_vtk_write_cart_rho_f64(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_rho_f64.vtk"
    vtkdata_cart.export_to_vtk(filename=file, time=0, names="rho", dtype="float64")
    expected = tmpdir / "test_cart_rho_f64_t0000.vtk"
    assert expected.is_file()


def test_vtk_write_invalid_dtype(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_rho_invalid_dtype.vtk"
    with pytest.raises(ValueError):
        vtkdata_cart.export_to_vtk(filename=file, time=0, names="rho", dtype="int32")


def test_vtk_write_multiple(vtkdata_cyl, tmpdir):
    file = tmpdir / "test_cyl_multiple.vtk"
    time = np.arange(5)
    vtkdata_cyl.export_to_vtk(
        filename=file, time=time, names=["rho", "v1", "v2"], bg_names=["rho0", "T0"]
    )
    for t in time:
        expected = tmpdir / f"test_cyl_multiple_t000{t}.vtk"
        assert expected.is_file()


def test_vtk_write_derived(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_derived.vtk"
    vtkdata_cart.export_to_vtk(
        filename=file, time=0, names=["rho", "v1", "B2", "B3"], bg_names=["rho0"]
    )
    expected = tmpdir / "test_cart_derived_t0000.vtk"
    assert expected.is_file()


def test_vtk_write_derived_not_present(vtkdata_cyl, tmpdir):
    file = tmpdir / "test_cyl_derived_not_present.vtk"
    with pytest.raises(ValueError):
        vtkdata_cyl.export_to_vtk(
            filename=file, time=0, names=["rho", "B2"], bg_names=["v02"]
        )


def test_vtk_write_invalid_ef(vtkdata_cart, tmpdir):
    file = tmpdir / "test_cart_invalid_ef.vtk"
    with pytest.raises(ValueError):
        vtkdata_cart.export_to_vtk(filename=file, time=0, names=["rho1"])
