"""AMICI-generated module for model Smith_BMCSystBiol2013"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Smith_BMCSystBiol2013` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Smith_BMCSystBiol2013 = amici._module_from_path(
    "Smith_BMCSystBiol2013.Smith_BMCSystBiol2013", Path(__file__).parent / "Smith_BMCSystBiol2013.py"
)
for var in dir(Smith_BMCSystBiol2013):
    if not var.startswith("_"):
        globals()[var] = getattr(Smith_BMCSystBiol2013, var)
get_model = Smith_BMCSystBiol2013.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Smith_BMCSystBiol2013._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Smith_BMCSystBiol2013._model_module = sys.modules[__name__]

__version__ = "0.1.0"
