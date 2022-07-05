
import os

from xanthos import Xanthos


def run(ini):

    # instantiate model
    xth = Xanthos(ini)

    # run intended configuration
    xth.execute()

    return xth


if __name__ == "__main__":

    # full path to parameterized config file
    ini = os.path.join('/workflow/runoff-watch-setup/xanthos76/pm_abcd_mrtm76.ini')

    # run the model
    xth = run(ini)

    del xth