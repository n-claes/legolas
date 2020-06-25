import pylbo

FP_LIMIT = 1e-12
output = (pylbo.LEGOLAS_DIR / 'tests/regression_tests').resolve()


def get_filepaths(basename):
    filename_dat = '{}.dat'.format(basename)
    filename_log = '{}.log'.format(basename)
    datfile = (output / filename_dat).resolve()
    logfile = (output / filename_log).resolve()
    return datfile, logfile


def get_answer_filepaths(basename):
    datfile, logfile = get_filepaths(basename)
    answer_datfile = (output / 'answers/answer_{}'.format(datfile.name)).resolve()
    answer_logfile = (output / 'answers/answer_{}'.format(logfile.name)).resolve()
    return answer_datfile, answer_logfile
