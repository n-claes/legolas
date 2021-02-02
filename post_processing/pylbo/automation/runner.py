import os
import subprocess
import signal
import multiprocessing
import psutil
import tqdm
from pathlib import Path

from pylbo.utilities.logger import pylboLogger
from pylbo.utilities.toolbox import transform_to_list


def _validate_executable(executable):
    """
    Validates the given executable, then returns it. If the argument passed is
    None, defaults to the executable in the legolas home directory.

    Parameters
    ----------
    executable : str, ~os.PathLike
        The path to the legolas executable to use.

    Raises
    ------
    NotADirectoryError
        If the directory containing the executable was not found.
    FileNotFoundError
        If the executable was not found.

    Returns
    -------
    executable : ~os.PathLike
        The (resolved) path to the executable to use.

    """
    if executable is None:
        exec_dir = Path(__file__).parents[3]
        executable = (exec_dir / "legolas").resolve()
    executable = Path(executable).resolve()
    if not executable.parent.is_dir():
        raise NotADirectoryError(
            f"Directory containing the executable was not found: {executable.parent}"
        )
    if not executable.is_file():
        raise FileNotFoundError(f"Executable was not found: {executable}")
    return executable


def _validate_nb_cpus(cpus):
    """
    Validates the number of cpus passed to the multiprocessing pool.
    Defaults to the maximum available number if exceeded.

    Parameters
    ----------
    cpus : int
        The number of cpus to use.

    Returns
    -------
    cpus : int
        The number of cpus to use, limited to the maximum number available.

    """
    cpus_available = multiprocessing.cpu_count()
    if cpus > cpus_available:
        pylboLogger.warning(
            f"Requested more than the available number of cpus ({cpus}). "
            f"Setting nb_cpus to maximum available ({cpus_available})."
        )
        cpus = cpus_available
    return cpus


def _validate_parfiles(files):
    """
    Validates a list of parfiles.

    Parameters
    ----------
    files : (list of) str, (list of) ~os.PathLike
        Paths to the parfiles.

    Raises
    ------
    FileNotFoundError
        If one of the parfiles is not found.

    Returns
    -------
    files_list : list
        A list of resolved filepaths for the parfiles.

    """
    files_list = transform_to_list(files)
    files_list = [Path(file).resolve() for file in files_list]
    for file in files_list:
        if not file.is_file():
            raise FileNotFoundError(f"Parfile was not found: {file}")
    return files_list


class LegolasRunner:
    """
    Handles running legolas.

    Parameters
    ----------
    parfiles : list, numpy.ndarray
        A list or array containing the names or paths to the parfiles.
    remove_parfiles : bool
        If `True`, removes the parfiles after running Legolas. This will also remove
        the containing folder if it turns out to be empty after the parfiles are
        removed. If there are other files still in the folder it remains untouched.
    nb_cpus : int
        The number of CPUs to use when running Legolas. If equal to 1 then
        parallelisation is disabled. Defaults to the maximum number of CPUs available
        if a number larger than the available number is specified.
    executable : str, ~os.PathLike
        The path to the legolas executable. If not specified, defaults to the
        standard one in the legolas home directory.
    """

    def __init__(self, parfiles, remove_parfiles, nb_cpus, executable):
        self.parfiles = _validate_parfiles(parfiles)
        self.parfile_dir = self.parfiles[0].parent
        self.executable = _validate_executable(executable)
        self.nb_cpus = _validate_nb_cpus(nb_cpus)
        self.remove_parfiles = remove_parfiles

        pylboLogger.info(f"initialising runner, using executable {self.executable}")

    @staticmethod
    def _init_worker():
        """
        Worker initialisation for the multiprocessing module.
        """
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    def _activate_worker(self, parfile):
        """
        Worker activation for the multiprocessing module.
        Calls the legolas executable as a subprocess with the parfile as argument.

        Parameters
        ----------
        parfile : str, ~os.PathLike
            The path to the parfile

        Returns
        -------
        call : :func:`subprocess.call`
            A call to a subprocess to run legolas.
        """
        cmd = [f"./{self.executable.stem}", "-i", str(parfile)]
        return subprocess.call(cmd)

    @staticmethod
    def _terminate_workers():
        """
        Terminates the multiprocessing workers after a forced interruption.
        Simply giving an interrupt terminates only the Python processes, but still
        keeps the Legolas calls running since those are subprocesses. This method first
        terminates all child processes (legolas), then the parents (workers).
        """
        pylboLogger.error("interrupting processes...")
        for process in multiprocessing.active_children():
            pid = process.pid
            pylboLogger.error(f"terminating PID: {pid} -- {process.name}")
            parent = psutil.Process(pid)
            children = parent.children(recursive=True)
            for child in children:
                pylboLogger.error(
                    f"terminating child process {child.pid} -- {child.name}"
                )
                child.kill()
            gone, alive = psutil.wait_procs(children, timeout=2)
            for killed_proc in gone:
                pylboLogger.error(f"{str(killed_proc)}")
            parent.kill()
            parent.wait(timeout=2)
        pylboLogger.critical("all Legolas processes terminated.")
        exit(1)

    def execute(self):
        """
        Executes the legolas executables and initialises the multiprocessing
        pool if requested.
        """

        def update_pbar(*args):
            pbar.update()

        # original working directory
        owd = Path.cwd()
        # change to parent directory of executable
        os.chdir(self.executable.parent)
        # no parallelisation if there is only one parfile
        if len(self.parfiles) == 1:
            try:
                pylboLogger.info("running legolas...")
                self._activate_worker(*self.parfiles)
            except KeyboardInterrupt:
                self._terminate_workers()
                exit(1)
        else:
            # initialise progressbar and multiprocessing pool
            pbar = tqdm.tqdm(total=len(self.parfiles), unit="")
            pbar.set_description(f"running legolas [{self.nb_cpus} CPUS]")
            pool = multiprocessing.Pool(
                processes=self.nb_cpus, initializer=self._init_worker
            )
            try:
                for parfile in self.parfiles:
                    pool.apply_async(
                        self._activate_worker, args=(parfile,), callback=update_pbar
                    )
                pool.close()
                pool.join()
                pbar.close()
            except KeyboardInterrupt:
                pbar.set_description("INTERRUPTED")
                pbar.update(len(self.parfiles))
                pbar.close()
                self._terminate_workers()
                pool.terminate()
                pool.join()
                exit(1)
        pylboLogger.info("all runs completed")

        # change back to the original directory
        os.chdir(owd)
        # if requested, remove parfiles
        if self.remove_parfiles:
            for file in self.parfiles:
                os.remove(file)
            pylboLogger.info("parfiles removed.")
            # if directory is empty, also remove it
            try:
                Path.rmdir(self.parfile_dir)
                pylboLogger.info(
                    f"parfile containing folder '{self.parfile_dir}' "
                    f"was empty and is also removed."
                )
            except OSError:
                pass
