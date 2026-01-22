"""AMICI-generated module for model Borghans_BiophysChem1997"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Borghans_BiophysChem1997` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Borghans_BiophysChem1997 = amici._module_from_path(
    "Borghans_BiophysChem1997.Borghans_BiophysChem1997", Path(__file__).parent / "Borghans_BiophysChem1997.py"
)
for var in dir(Borghans_BiophysChem1997):
    if not var.startswith("_"):
        globals()[var] = getattr(Borghans_BiophysChem1997, var)
get_model = Borghans_BiophysChem1997.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Borghans_BiophysChem1997._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Borghans_BiophysChem1997._model_module = sys.modules[__name__]

__version__ = "0.1.0"
