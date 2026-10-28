"""Publish explicitly allowlisted assets from a version tag using an Actions token.

No credentials are stored or printed. Existing published releases are not modified.
"""
import hashlib
import json
import os
from pathlib import Path
import runpy
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen


def main():
    root = Path(__file__).resolve().parent.parent
    version = runpy.run_path(str(root / "idgrm/_version.py"))["__version__"]
    tag = os.environ["GITHUB_REF_NAME"]
    repository = os.environ["GITHUB_REPOSITORY"]
    if tag != f"v{version}" or repository != "lihongxiang97/iDGRM":
        raise ValueError("Release tag, version or repository does not match the intended publication")
    token = os.environ["GITHUB_TOKEN"]
    base = f"https://api.github.com/repos/{repository}"

    def request(url, method="GET", data=None, content_type="application/json"):
        headers = {"Authorization": "Bearer " + token, "Accept": "application/vnd.github+json",
            "User-Agent": "iDGRM-release", "Content-Type": content_type}
        encoded = json.dumps(data).encode() if data is not None and content_type == "application/json" else data
        with urlopen(Request(url, data=encoded, method=method, headers=headers), timeout=60) as response:
            return json.load(response)

    assets = [root / "dist" / name for name in (
        f"idgrm-{version}-py3-none-any.whl", f"iDGRM-{version}-source.zip",
        f"iDGRM-{version}-source.tar.gz", f"SHA256SUMS_{version}.json")]
    checksums = json.loads(assets[-1].read_text(encoding="utf-8"))
    for path in assets[:-1]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != checksums[path.name]:
            raise ValueError(f"Asset checksum mismatch: {path.name}")
    try:
        release = request(base + "/releases/tags/" + quote(tag))
    except HTTPError as error:
        if error.code != 404:
            raise
        release = request(base + "/releases", "POST", dict(tag_name=tag, name=f"iDGRM {version}",
            body=(root / f"docs/RELEASE_{version}_CN.md").read_text(encoding="utf-8"),
            draft=True, prerelease=False))
    if not release["draft"]:
        raise ValueError("Release is already published; refusing to overwrite")
    if release["assets"]:
        raise ValueError("Draft already contains assets; refusing implicit replacement")
    upload = release["upload_url"].split("{", 1)[0]
    for path in assets:
        result = request(upload + "?name=" + quote(path.name), "POST", path.read_bytes(), "application/octet-stream")
        print("Uploaded", result["name"], result["size"])
    published = request(base + "/releases/" + str(release["id"]), "PATCH", {"draft": False})
    print("Published", published["html_url"])


if __name__ == "__main__":
    main()
