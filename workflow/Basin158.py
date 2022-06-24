
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
    ini = os.path.join('/project/hli/gabeshu/Guta_Working/BasinsFile/xanthos158/pm_abcd_mrtm158.ini')

    # run the model
    xth = run(ini)

    del xth