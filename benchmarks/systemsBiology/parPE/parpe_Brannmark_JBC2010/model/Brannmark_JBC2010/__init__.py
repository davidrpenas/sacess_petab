"""AMICI-generated module for model Brannmark_JBC2010"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Brannmark_JBC2010` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Brannmark_JBC2010 = amici._module_from_path(
    "Brannmark_JBC2010.Brannmark_JBC2010", Path(__file__).parent / "Brannmark_JBC2010.py"
)
for var in dir(Brannmark_JBC2010):
    if not var.startswith("_"):
        globals()[var] = getattr(Brannmark_JBC2010, var)
get_model = Brannmark_JBC2010.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Brannmark_JBC2010._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Brannmark_JBC2010._model_module = sys.modules[__name__]

__version__ = "0.1.0"
