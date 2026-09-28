import numpy as np
import pytest

from spectro.core.io import (ParseOptions, export_csv, export_excel, load_file, write_jcamp,
                             write_spc)
from spectro.core.io.jcamp import decode_asdf_line, parse_jcamp
from spectro.core.spectrum import Spectrum

X = np.arange(200, 401, 1.0)


def gauss(c=270, w=20):
    return np.exp(-(((X - c) / w) ** 2))


def test_single_column_with_title_and_descending_axis(tmp_path):
    f = tmp_path / "a.txt"
    f.write_text('"Sample A"\n"Wavelength nm.","Abs."\n' + "\n".join(
        f"{a:.1f},{v:.5f}" for a, v in zip(X[::-1], gauss()[::-1])))
    r = load_file(f)
    s = r.spectra[0]
    assert r.layout == "columns" and s.name == "Sample A"
    assert s.x[0] == 200 and np.allclose(s.y, gauss(), atol=1e-5)


def test_semicolon_decimal_comma_multicolumn_and_metadata(tmp_path):
    f = tmp_path / "b.csv"
    f.write_text("Operator;Ali\nWavelength;S1;S2\n" + "\n".join(
        f"{a:.1f};0.500;0.250".replace(".", ",") for a in X))
    r = load_file(f)
    assert [s.name for s in r.spectra] == ["S1", "S2"]
    assert r.spectra[1].y[0] == 0.25 and r.metadata["Operator"] == "Ali"


def test_xy_pairs_layout(tmp_path):
    lines = ["Std1,,Std2,,", "Wavelength (nm),Abs,Wavelength (nm),Abs,"]
    lines += [f"{a},0.1,{a},0.2," for a in X[::-1]]
    f = tmp_path / "c.csv"
    f.write_text("\n".join(lines))
    r = load_file(f)
    assert r.layout == "xy_pairs" and [s.name for s in r.spectra] == ["Std1", "Std2"]


def test_rows_layout_tab(tmp_path):
    f = tmp_path / "d.tsv"
    f.write_text("Sample\t" + "\t".join(map(str, X)) + "\nmix1\t" + "\t".join(["0.3"] * X.size)
                 + "\nmix2\t" + "\t".join(["0.4"] * X.size))
    r = load_file(f)
    assert r.layout == "rows" and [s.name for s in r.spectra] == ["mix1", "mix2"]


def test_whitespace_with_comment_header(tmp_path):
    f = tmp_path / "e.asc"
    f.write_text("# exported\n" + "\n".join(f"{a}   0.7" for a in X))
    r = load_file(f)
    assert r.spectra[0].name == "e" and r.spectra[0].x.size == X.size


def test_utf16_file(tmp_path):
    f = tmp_path / "u.txt"
    f.write_bytes("nm\tAbs\n".encode("utf-16") + "".join(
        f"{a}\t{v:.4f}\n" for a, v in zip(X, gauss())).encode("utf-16-le"))
    r = load_file(f)
    assert r.spectra[0].x.size == X.size


def test_forced_layout(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("\n".join(f"{a},{v}" for a, v in zip(X, gauss())))
    with pytest.raises(ValueError):
        load_file(f, ParseOptions(layout="rows"))


def test_roundtrips(tmp_path):
    S = [Spectrum(X, np.sin(X / 30) + 2, name="s1"), Spectrum(X, np.cos(X / 30) + 2, name="s,2")]
    write_jcamp(S[0], tmp_path / "a.jdx")
    assert np.allclose(load_file(tmp_path / "a.jdx").spectra[0].y, S[0].y, atol=1e-6)
    write_spc(S, tmp_path / "a.spc")
    r = load_file(tmp_path / "a.spc")
    assert len(r.spectra) == 2 and np.allclose(r.spectra[1].y, S[1].y, atol=1e-5)
    export_excel(S, tmp_path / "a.xlsx")
    r = load_file(tmp_path / "a.xlsx")
    assert [s.name for s in r.spectra] == ["s1", "s,2"]
    export_csv(S, tmp_path / "b.csv")
    assert [s.name for s in load_file(tmp_path / "b.csv").spectra] == ["s1", "s,2"]


def test_jcamp_compressed_forms():
    assert decode_asdf_line("8192@jJ%kJ") == [8192, 0, -1, 0, 0, -2, -1]
    assert decode_asdf_line("1@JT") == [1, 0, 1, 2]
    text = """##TITLE=test
##JCAMP-DX=5.01
##XUNITS=NANOMETERS
##YUNITS=ABSORBANCE
##XFACTOR=1
##YFACTOR=0.001
##FIRSTX=200
##LASTX=207
##DELTAX=1
##NPOINTS=8
##XYDATA=(X++(Y..Y))
200@A0J0J0J0
204D0J0J0J0
##END="""
    s = parse_jcamp(text).spectra[0]
    assert np.allclose(s.x, np.arange(200, 208))
    assert np.allclose(s.y, [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])


def test_tab_delimited_rows_with_spaces_in_names(tmp_path):
    f = tmp_path / "v.tsv"
    f.write_text("Sample\t" + "\t".join(map(str, X)) + "\nV1 PAR 9 CAF 5\t"
                 + "\t".join(["0.3"] * X.size))
    r = load_file(f)
    assert r.spectra[0].name == "V1 PAR 9 CAF 5" and r.spectra[0].x.size == X.size
