"""AMICI-generated module for model Raia_CancerResearch2011"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Raia_CancerResearch2011` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Raia_CancerResearch2011 = amici._module_from_path(
    "Raia_CancerResearch2011.Raia_CancerResearch2011", Path(__file__).parent / "Raia_CancerResearch2011.py"
)
for var in dir(Raia_CancerResearch2011):
    if not var.startswith("_"):
        globals()[var] = getattr(Raia_CancerResearch2011, var)
get_model = Raia_CancerResearch2011.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Raia_CancerResearch2011._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Raia_CancerResearch2011._model_module = sys.modules[__name__]

__version__ = "0.1.0"
