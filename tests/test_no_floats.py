"""Structural enforcement of the loud-exactness contract (spec D3, §4.1):
no float/complex literals and no float() calls anywhere under src/."""
import ast
import pathlib

SRC = pathlib.Path(__file__).resolve().parent.parent / "src" / "quiverlab"


def _violations(path: pathlib.Path) -> list[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)):
            out.append(f"{path.name}:{node.lineno}: literal {node.value!r}")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "float"
        ):
            out.append(f"{path.name}:{node.lineno}: float() call")
    return out


def test_package_importable():
    import quiverlab
    assert isinstance(quiverlab.__version__, str)


def test_no_float_literals_or_calls_in_src():
    assert SRC.is_dir(), "src/quiverlab missing"
    bad = [v for f in SRC.rglob("*.py") for v in _violations(f)]
    assert bad == [], "floats are banned in quiverlab core:\n" + "\n".join(bad)
