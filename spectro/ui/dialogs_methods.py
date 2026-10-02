"""Quantitative method dialogs (compatibility module).

The dialogs live in :mod:`spectro.ui.methods`, one module per method family;
everything is re-exported here so existing imports keep working."""

from PySide6.QtWidgets import QFileDialog, QInputDialog  # noqa: F401

from spectro.storage.recalc import determine, method_from_definition  # noqa: F401
from spectro.ui.methods.chemometrics import ChemometricsDialog  # noqa: F401
from spectro.ui.methods.common import (CAL_ROLES, SURFACE_RING, TEST_ROLES,  # noqa: F401
                                       TRAIN_ROLES, CompoundChecks, add_grouped, check_ids,
                                       checklists, list_box, results_rows, scatter, scroll,
                                       start_editing, store_method, summary_text, xy_plot)
from spectro.ui.methods.equations import EquationsDialog  # noqa: F401
from spectro.ui.methods.progressive import (ProgressiveDialog, _cell, _num_cell,  # noqa: F401
                                            pure_standards)
from spectro.ui.methods.saved import ResultEditDialog, SavedDialog  # noqa: F401
from spectro.ui.methods.special import SpecialDialog, wl_spin  # noqa: F401
from spectro.ui.methods.univariate import RobustnessDialog, UnivariateDialog  # noqa: F401
