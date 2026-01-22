"""AMICI-generated module for model Zheng_PNAS2012"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Zheng_PNAS2012` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Zheng_PNAS2012 = amici._module_from_path(
    "Zheng_PNAS2012.Zheng_PNAS2012", Path(__file__).parent / "Zheng_PNAS2012.py"
)
for var in dir(Zheng_PNAS2012):
    if not var.startswith("_"):
        globals()[var] = getattr(Zheng_PNAS2012, var)
get_model = Zheng_PNAS2012.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Zheng_PNAS2012._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Zheng_PNAS2012._model_module = sys.modules[__name__]

__version__ = "0.1.0"
