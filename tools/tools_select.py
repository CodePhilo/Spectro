"""Select spectra in the main window's project tree (screenshot tools)."""

from PySide6.QtCore import Qt


def select_ids(win, ids) -> None:
    win.tree.clearSelection()
    stack = [win.tree.invisibleRootItem()]
    while stack:
        node = stack.pop()
        for i in range(node.childCount()):
            ch = node.child(i)
            if ch.data(0, Qt.UserRole) in ids:
                ch.setSelected(True)
            ch.setExpanded(True)
            stack.append(ch)
