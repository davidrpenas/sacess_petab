"""AMICI-generated module for model Beer_MolBioSystems2014"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Beer_MolBioSystems2014` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Beer_MolBioSystems2014 = amici._module_from_path(
    "Beer_MolBioSystems2014.Beer_MolBioSystems2014", Path(__file__).parent / "Beer_MolBioSystems2014.py"
)
for var in dir(Beer_MolBioSystems2014):
    if not var.startswith("_"):
        globals()[var] = getattr(Beer_MolBioSystems2014, var)
get_model = Beer_MolBioSystems2014.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Beer_MolBioSystems2014._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Beer_MolBioSystems2014._model_module = sys.modules[__name__]

__version__ = "0.1.0"
