"""AMICI-generated module for model Sneyd_PNAS2002"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Sneyd_PNAS2002` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Sneyd_PNAS2002 = amici._module_from_path(
    "Sneyd_PNAS2002.Sneyd_PNAS2002", Path(__file__).parent / "Sneyd_PNAS2002.py"
)
for var in dir(Sneyd_PNAS2002):
    if not var.startswith("_"):
        globals()[var] = getattr(Sneyd_PNAS2002, var)
get_model = Sneyd_PNAS2002.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Sneyd_PNAS2002._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Sneyd_PNAS2002._model_module = sys.modules[__name__]

__version__ = "0.1.0"
