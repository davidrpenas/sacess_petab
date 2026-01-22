"""AMICI-generated module for model Chen_MSB2009"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Chen_MSB2009` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Chen_MSB2009 = amici._module_from_path(
    "Chen_MSB2009.Chen_MSB2009", Path(__file__).parent / "Chen_MSB2009.py"
)
for var in dir(Chen_MSB2009):
    if not var.startswith("_"):
        globals()[var] = getattr(Chen_MSB2009, var)
get_model = Chen_MSB2009.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Chen_MSB2009._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Chen_MSB2009._model_module = sys.modules[__name__]

__version__ = "0.1.0"
