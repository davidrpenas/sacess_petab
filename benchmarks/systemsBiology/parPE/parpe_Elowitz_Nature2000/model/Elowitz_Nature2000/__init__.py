"""AMICI-generated module for model Elowitz_Nature2000"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Elowitz_Nature2000` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Elowitz_Nature2000 = amici._module_from_path(
    "Elowitz_Nature2000.Elowitz_Nature2000", Path(__file__).parent / "Elowitz_Nature2000.py"
)
for var in dir(Elowitz_Nature2000):
    if not var.startswith("_"):
        globals()[var] = getattr(Elowitz_Nature2000, var)
get_model = Elowitz_Nature2000.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Elowitz_Nature2000._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Elowitz_Nature2000._model_module = sys.modules[__name__]

__version__ = "0.1.0"
