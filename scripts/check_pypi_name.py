"""Advisory: is a name free on PyPI? 404 on the JSON API == available.
(The /project/<name>/ page is Cloudflare-challenged and unreliable via curl.)"""
import urllib.error
import urllib.request


def available(name: str):
    try:
        urllib.request.urlopen(f"https://pypi.org/pypi/{name}/json", timeout=10)
        return False                       # 200 -> taken
    except urllib.error.HTTPError as e:
        return True if e.code == 404 else None
    except Exception:                      # noqa: BLE001  (offline -> unknown)
        return None


if __name__ == "__main__":
    for n in ("quiverlab", "quiver-lab"):
        a = available(n)
        print(f"{n}: {'AVAILABLE' if a else 'TAKEN' if a is False else 'UNKNOWN (offline?)'}")
