"""AMICI-generated module for model Bruno_JExpBot2016"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Bruno_JExpBot2016` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Bruno_JExpBot2016 = amici._module_from_path(
    "Bruno_JExpBot2016.Bruno_JExpBot2016", Path(__file__).parent / "Bruno_JExpBot2016.py"
)
for var in dir(Bruno_JExpBot2016):
    if not var.startswith("_"):
        globals()[var] = getattr(Bruno_JExpBot2016, var)
get_model = Bruno_JExpBot2016.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Bruno_JExpBot2016._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Bruno_JExpBot2016._model_module = sys.modules[__name__]

__version__ = "0.1.0"
