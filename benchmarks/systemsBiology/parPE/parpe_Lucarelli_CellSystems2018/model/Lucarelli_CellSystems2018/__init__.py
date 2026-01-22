"""AMICI-generated module for model Lucarelli_CellSystems2018"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Lucarelli_CellSystems2018` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Lucarelli_CellSystems2018 = amici._module_from_path(
    "Lucarelli_CellSystems2018.Lucarelli_CellSystems2018", Path(__file__).parent / "Lucarelli_CellSystems2018.py"
)
for var in dir(Lucarelli_CellSystems2018):
    if not var.startswith("_"):
        globals()[var] = getattr(Lucarelli_CellSystems2018, var)
get_model = Lucarelli_CellSystems2018.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Lucarelli_CellSystems2018._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Lucarelli_CellSystems2018._model_module = sys.modules[__name__]

__version__ = "0.1.0"
