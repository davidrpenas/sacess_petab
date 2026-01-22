"""AMICI-generated module for model Armistead_CellDeathDis2024"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Armistead_CellDeathDis2024` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Armistead_CellDeathDis2024 = amici._module_from_path(
    "Armistead_CellDeathDis2024.Armistead_CellDeathDis2024", Path(__file__).parent / "Armistead_CellDeathDis2024.py"
)
for var in dir(Armistead_CellDeathDis2024):
    if not var.startswith("_"):
        globals()[var] = getattr(Armistead_CellDeathDis2024, var)
get_model = Armistead_CellDeathDis2024.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Armistead_CellDeathDis2024._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Armistead_CellDeathDis2024._model_module = sys.modules[__name__]

__version__ = "0.1.0"
