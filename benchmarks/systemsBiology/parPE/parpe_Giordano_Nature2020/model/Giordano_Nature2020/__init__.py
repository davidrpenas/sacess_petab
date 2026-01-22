"""AMICI-generated module for model Giordano_Nature2020"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Giordano_Nature2020` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Giordano_Nature2020 = amici._module_from_path(
    "Giordano_Nature2020.Giordano_Nature2020", Path(__file__).parent / "Giordano_Nature2020.py"
)
for var in dir(Giordano_Nature2020):
    if not var.startswith("_"):
        globals()[var] = getattr(Giordano_Nature2020, var)
get_model = Giordano_Nature2020.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Giordano_Nature2020._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Giordano_Nature2020._model_module = sys.modules[__name__]

__version__ = "0.1.0"
