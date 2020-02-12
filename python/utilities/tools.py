import sys

def custom_enumerate(iterable, start=0, step=1):
    for itr in iterable:
        yield start, itr
        start += step


def progressbar(current, total, text=''):
    bar_length = 40
    progress = float(current) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(bar_length * progress))
    text = "\r[{}] {:.0f}% {}".format("#" * block + "-" * (bar_length - block), round(progress * 100, 0), text)
    sys.stdout.write(text)
    sys.stdout.flush()

