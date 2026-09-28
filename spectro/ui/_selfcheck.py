"""Import every UI module and core dependency (used to verify frozen builds)."""


def run() -> str:
    import olefile
    import matplotlib.backends.backend_agg
    import openpyxl
    import sklearn.cross_decomposition
    import sklearn.neural_network
    import sklearn.svm
    import xlrd

    from spectro.demo import data_dir
    from spectro.ui import dialogs_data, dialogs_methods, dialogs_tools, figure, report

    if not (data_dir() / "ternary" / "excedrin_tablets.spc").exists():
        raise RuntimeError("demo data missing from the build")

    modules = (matplotlib.backends.backend_agg, olefile, openpyxl, sklearn.cross_decomposition, sklearn.neural_network,
               sklearn.svm, xlrd, dialogs_data, dialogs_methods, dialogs_tools, figure,
               report)
    return f"selfcheck ok ({len(modules)} modules)"
