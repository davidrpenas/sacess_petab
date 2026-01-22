"""AMICI-generated module for model Schwen_PONE2014"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Schwen_PONE2014` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Schwen_PONE2014 = amici._module_from_path(
    "Schwen_PONE2014.Schwen_PONE2014", Path(__file__).parent / "Schwen_PONE2014.py"
)
for var in dir(Schwen_PONE2014):
    if not var.startswith("_"):
        globals()[var] = getattr(Schwen_PONE2014, var)
get_model = Schwen_PONE2014.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Schwen_PONE2014._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Schwen_PONE2014._model_module = sys.modules[__name__]

__version__ = "0.1.0"
