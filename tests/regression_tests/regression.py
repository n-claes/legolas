import pytest
import pylbo
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
from pathlib import Path

testlog = logging.getLogger("test_logger")

only_for_baseline_generation = pytest.mark.skipif(
    condition="not config.getoption('generate_baseline')",
    reason="'--generate' option not passed",
)


class TestCase:
    SAVEFIG_KWARGS = {"dpi": 200, "transparent": True}
    RMS_TOLERANCE = 2

    gridpoints = 51
    logging_level = 1
    show_results = False

    @property
    def name(self):
        raise NotImplementedError()

    @property
    def filename(self):
        raise NotImplementedError()

    @property
    def parameters(self):
        raise NotImplementedError()

    @property
    def equilibrium(self):
        raise NotImplementedError()

    @property
    def geometry(self):
        raise NotImplementedError()

    @property
    def eigenfunction_settings(self):
        return {"write_eigenfunctions": False}

    @property
    def physics_settings(self):
        return {}

    @property
    def eigenvalues_are_real(self):
        return False

    def get_spectrum_image_filenames(self, limits):
        xlim = limits["xlim"]
        ylim = limits["ylim"]
        figname_test = f"{self.filename}_Re{xlim[0]}-{xlim[1]}_Im{ylim[0]}-{ylim[1]}"
        figname_base = f"{figname_test}-baseline"
        return (
            self._spectradir / f"{figname_test}.png",
            self._spectradir / f"{figname_base}.png",
        )

    def get_eigenfunction_image_filenames(self, eigenvalue):
        figname_test = f"{self.filename}_efs_w_{eigenvalue:.8f}"
        figname_base = f"{figname_test}-baseline"
        return (
            self._eigfuncdir / f"{figname_test}.png",
            self._eigfuncdir / f"{figname_base}.png",
        )

    @pytest.fixture(scope="class")
    def file_base(self, baselinedir):
        return baselinedir / f"BASE_{self.filename}.dat"

    @pytest.fixture(scope="class")
    def file_test(self, datfiledir):
        return datfiledir / f"{self.filename}.dat"


class RegressionTest(TestCase):
    def setup(self, outputdir):
        _setup = {
            "geometry": self.geometry,
            "x_start": getattr(self, "x_start", 0),
            "x_end": getattr(self, "x_end", 1),
            "gridpoints": self.gridpoints,
            "parameters": self.parameters,
            "equilibrium_type": self.equilibrium,
            "logging_level": self.logging_level,
            "show_results": self.show_results,
            "basename_datfile": self.filename,
            "output_folder": str(outputdir),
        }
        _setup.update(self.eigenfunction_settings)
        _setup.update(self.physics_settings)
        return _setup

    @only_for_baseline_generation
    def test_generate_baseline(self, pytestconfig, file_base):
        if file_base.is_file():
            testlog.info(f"baseline file '{file_base.name}' is already present!")
            capturemanager = pytestconfig.pluginmanager.getplugin("capturemanager")
            capturemanager.suspend_global_capture(in_=True)
            override = input("Regenerate (replace) this file? ").lower() in ("yes", "y")
            capturemanager.resume_global_capture()
            if not override:
                pytest.skip("using existing file")

        setup = self.setup(outputdir=file_base.parent)
        setup.update({"basename_datfile": file_base.stem})
        self.generate_test_dataset(setup)

    @pytest.fixture(scope="class")
    def ds_base(self, file_base):
        return pylbo.load(file_base)

    @pytest.fixture(scope="class")
    def ds_test(self, file_test, datfiledir):
        setup = self.setup(datfiledir)
        self.generate_test_dataset(setup)
        return pylbo.load(file_test)

    def generate_test_dataset(self, setup):
        testlog.info(f"generating dataset: {setup['basename_datfile']}.dat")
        parfile = pylbo.generate_parfiles(
            parfile_dict=setup,
            basename=setup["basename_datfile"],
            output_dir=self._datfiledir,
            subdir=False,
        )
        pylbo.run_legolas(parfile)

    def generate_spectrum_images(self, limits, ds_test, ds_base):
        p_test = pylbo.plot_spectrum(ds_test)
        p_base = pylbo.plot_spectrum(ds_base)
        figname_test, figname_base = self.get_spectrum_image_filenames(limits)
        xlim = limits["xlim"]
        ylim = limits["ylim"]
        for pp, name in [(p_test, figname_test), (p_base, figname_base)]:
            pp.ax.set_xlim(xlim)
            pp.ax.set_ylim(ylim)
            pp.ax.set_title(self.name)
            pp.fig.savefig(name, **self.SAVEFIG_KWARGS)
            plt.close(pp.fig)
        return (figname_test, figname_base)

    def generate_eigenfunction_images(self, eigenvalue, ds_test, ds_base):
        fig_test, ax_test = plt.subplots(3, 3, figsize=(10, 10), sharex="all")
        fig_base, ax_base = plt.subplots(3, 3, figsize=(10, 10), sharex="all")
        figname_test, figname_base = self.get_eigenfunction_image_filenames(eigenvalue)
        for ds, ax in [(ds_test, ax_test), (ds_base, ax_base)]:
            (efs,) = ds.get_eigenfunctions(eigenvalue)
            for panel, ef_name in zip(ax.flatten(), ds.ef_names):
                result = abs(efs[ef_name].real + efs[ef_name].imag)
                # small values
                result[np.where(result < 1e-10)] = 0
                panel.plot(ds.ef_grid, result, lw=3)
                panel.set_yticks([])
                panel.set_title(ef_name)
        for fig, name in [(fig_test, figname_test), (fig_base, figname_base)]:
            fig.suptitle(f"eigenvalue = {eigenvalue:.9f}")
            fig.tight_layout()
            fig.savefig(name, **self.SAVEFIG_KWARGS)
            plt.close(fig)
        return (figname_test, figname_base)

    def compare_test_images(self, image_test, image_baseline, tol):
        result = compare_images(str(image_baseline), str(image_test), tol=tol)
        if result is not None:
            pytest.fail(result, pytrace=False)
        # test succeeded if result = None, check if files are kept
        if result is None and not self._keep_files:
            Path(image_baseline).unlink()
            Path(image_test).unlink()

    def test_generate_ds(self, ds_test):
        assert ds_test is not None

    def test_file_base_exists(self, file_base):
        assert file_base.is_file()

    def test_ds_base_exists(self, ds_base):
        assert ds_base is not None

    def test_eigenvalue_types(self, ds_test):
        if self.eigenvalues_are_real:
            assert np.all(ds_test.eigenvalues.imag == pytest.approx(0))
        else:
            assert np.any(ds_test.eigenvalues.imag != pytest.approx(0))

    def test_geometry(self, ds_test, ds_base):
        assert self.geometry == ds_test.geometry == ds_base.geometry

    def test_resolution(self, ds_test, ds_base):
        assert ds_test.gridpoints == ds_base.gridpoints

    def run_spectrum_test(self, limits, ds_test, ds_base):
        image_test, image_baseline = self.generate_spectrum_images(
            limits, ds_test, ds_base
        )
        self.compare_test_images(
            image_test,
            image_baseline,
            tol=limits.get("RMS_TOLERANCE", self.RMS_TOLERANCE),
        )

    def run_eigenfunction_test(self, eigenfunction, ds_test, ds_base):
        eigenvalue = eigenfunction["eigenvalue"]
        image_test, image_baseline = self.generate_eigenfunction_images(
            eigenvalue, ds_test, ds_base
        )
        self.compare_test_images(
            image_test,
            image_baseline,
            tol=eigenfunction.get("RMS_TOLERANCE", self.RMS_TOLERANCE),
        )
