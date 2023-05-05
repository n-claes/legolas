import logging

import numpy as np
import pytest

testlog = logging.getLogger("test_logger_pylbo")

only_for_baseline_generation = pytest.mark.skipif(
    condition="not config.getoption('generate_baseline')",
    reason="'--generate' option not passed",
)


def use_existing_baseline(capturemanager, baseline):
    use_existing = False
    if baseline.is_file():
        testlog.info(f"baseline file '{baseline.name}' is already present!")
        capturemanager.suspend_global_capture(in_=True)
        use_existing = input("Regenerate this file? ").lower() not in ("yes", "y")
        capturemanager.resume_global_capture()
    return use_existing


class ModeVizTest:
    @property
    def filename(self):
        raise NotImplementedError

    def cbar_matches(self, view, mode_solution):
        actual = (view.cbar.vmin, view.cbar.vmax)
        expected = (mode_solution.min(), mode_solution.max())
        return np.allclose(actual, expected)

    @pytest.fixture(scope="session")
    def capturemanager(self, pytestconfig):
        return pytestconfig.pluginmanager.getplugin("capturemanager")

    @only_for_baseline_generation
    def test_generate_mode_solution_baseline(
        self, capturemanager, view, modebaselinedir
    ):
        file = modebaselinedir / self.filename
        if use_existing_baseline(capturemanager, file):
            pytest.skip("using existing file")
        testlog.info(f"generating dataset: {self.filename}")
        np.save(file, view.solutions)

    @pytest.fixture(scope="class")
    def mode_solution(self, modebaselinedir):
        return np.load(modebaselinedir / self.filename)

    def test_solution_shape(self, view, mode_solution):
        assert view.solutions.shape == mode_solution.shape

    def test_mode_solution(self, view, mode_solution):
        assert np.allclose(view.solutions, mode_solution)
