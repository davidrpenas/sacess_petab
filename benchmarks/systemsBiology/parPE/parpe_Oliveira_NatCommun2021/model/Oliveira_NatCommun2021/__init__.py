"""AMICI-generated module for model Oliveira_NatCommun2021"""

import sys
from pathlib import Path
import amici


# Ensure we are binary-compatible, see #556
if "0.31.0" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Oliveira_NatCommun2021` in {Path(__file__).parent}, "
        "generated with amici==0.31.0, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.31.0 or re-import the model with the amici "
        "version currently installed."
    )

Oliveira_NatCommun2021 = amici._module_from_path(
    "Oliveira_NatCommun2021.Oliveira_NatCommun2021", Path(__file__).parent / "Oliveira_NatCommun2021.py"
)
for var in dir(Oliveira_NatCommun2021):
    if not var.startswith("_"):
        globals()[var] = getattr(Oliveira_NatCommun2021, var)
get_model = Oliveira_NatCommun2021.getModel

try:
    # _self: this module; will be set during import
    #  via amici.import_model_module
    Oliveira_NatCommun2021._model_module = _self  # noqa: F821
except NameError:
    # when the model package is imported via `import`
    Oliveira_NatCommun2021._model_module = sys.modules[__name__]

__version__ = "0.1.0"
