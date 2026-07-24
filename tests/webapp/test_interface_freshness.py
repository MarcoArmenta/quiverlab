from webapp.tools.check_library_interface import check


def test_library_interface_has_not_drifted():
    problems = check()
    assert problems == [], "Library interface drift:\n" + "\n".join(problems)
