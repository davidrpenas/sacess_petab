"""AMICI-generated module for model Bertozzi_PNAS2020"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Bertozzi_PNAS2020` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Bertozzi_PNAS2020 = amici._module_from_path(
    "Bertozzi_PNAS2020.Bertozzi_PNAS2020", Path(__file__).parent / "Bertozzi_PNAS2020.py"
)
for var in dir(Bertozzi_PNAS2020):
    if not var.startswith("_"):
        globals()[var] = getattr(Bertozzi_PNAS2020, var)
get_model = Bertozzi_PNAS2020.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Bertozzi_PNAS2020._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Bertozzi_PNAS2020._model_module = sys.modules[__name__]

__version__ = "0.1.0"
