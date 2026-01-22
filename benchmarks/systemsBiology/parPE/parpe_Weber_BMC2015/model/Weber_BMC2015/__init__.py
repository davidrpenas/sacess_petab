"""AMICI-generated module for model Weber_BMC2015"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Weber_BMC2015` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Weber_BMC2015 = amici._module_from_path(
    "Weber_BMC2015.Weber_BMC2015", Path(__file__).parent / "Weber_BMC2015.py"
)
for var in dir(Weber_BMC2015):
    if not var.startswith("_"):
        globals()[var] = getattr(Weber_BMC2015, var)
get_model = Weber_BMC2015.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Weber_BMC2015._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Weber_BMC2015._model_module = sys.modules[__name__]

__version__ = "0.1.0"
