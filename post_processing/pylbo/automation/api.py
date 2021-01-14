from pylbo.automation.generator import ParfileGenerator


def generate_parfiles(parfile_dict, basename=None, output_dir=None):
    pfgen = ParfileGenerator(parfile_dict, basename, output_dir)
    pfgen.create_namelist_from_dict()
    parfiles = pfgen.generate_parfiles()
    return parfiles


def run_legolas(parfiles, remove_parfiles=False, np_cpus=None, executable=None):
    pass
