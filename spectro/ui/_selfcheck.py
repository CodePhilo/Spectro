"""Import every UI module and core dependency (used to verify frozen builds)."""


def run() -> str:
    import olefile
    import openpyxl
    import sklearn.cross_decomposition
    import sklearn.neural_network
    import sklearn.svm
    import xlrd

    from spectro.ui import dialogs_data, dialogs_methods, dialogs_tools, report

    modules = (olefile, openpyxl, sklearn.cross_decomposition, sklearn.neural_network,
               sklearn.svm, xlrd, dialogs_data, dialogs_methods, dialogs_tools, report)
    return f"selfcheck ok ({len(modules)} modules)"
