"""AMICI-generated module for model Crauste_CellSystems2017"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Crauste_CellSystems2017` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Crauste_CellSystems2017 = amici._module_from_path(
    "Crauste_CellSystems2017.Crauste_CellSystems2017", Path(__file__).parent / "Crauste_CellSystems2017.py"
)
for var in dir(Crauste_CellSystems2017):
    if not var.startswith("_"):
        globals()[var] = getattr(Crauste_CellSystems2017, var)
get_model = Crauste_CellSystems2017.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Crauste_CellSystems2017._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Crauste_CellSystems2017._model_module = sys.modules[__name__]

__version__ = "0.1.0"
