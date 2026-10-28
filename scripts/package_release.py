"""Build a source-only delivery archive from an explicit allowlist; never include credentials/data."""
import hashlib
import json
from pathlib import Path
import tarfile
import zipfile
import runpy

root = Path(__file__).resolve().parent.parent
version = runpy.run_path(str(root / "idgrm/_version.py"))["__version__"]
out = root / "dist"
out.mkdir(exist_ok=True)
paths = [root / "README.md", root / "pyproject.toml"]
paths += [p for p in (root / "idgrm").rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix in {".py", ".R", ".tsv"}]
paths += list((root / "examples").glob("*.tsv"))
paths += [root / "docs" / name for name in ("METHODOLOGY_CN.md", "OUTGROUP_GUIDE_CN.md", f"RELEASE_{version}_CN.md")]
paths += [root / "scripts" / name for name in ("verify_release.py", "package_release.py")]
paths += [root / "scripts/publish_release.py", root / ".github/workflows/release.yml"]
# Include the self-contained core release tests, not unrelated research simulation harnesses.
paths += [root / "tests" / name for name in ("test_io.py", "test_legacy.py", "test_pipeline.py", "test_public_scope.py", "test_count_inputs.py", "test_polarization.py", "test_classifier_invariants.py")]
paths = sorted(set(paths))
for suffix in ("zip", "tar.gz"):
    destination = out / f"iDGRM-{version}-source.{suffix}"
    if suffix == "zip":
        with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in paths:
                archive.write(path, f"iDGRM-{version}/" + path.relative_to(root).as_posix())
    else:
        with tarfile.open(destination, "w:gz") as archive:
            for path in paths:
                archive.add(path, arcname=f"iDGRM-{version}/" + path.relative_to(root).as_posix())
manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob(f"*{version}*")
    if p.is_file() and not p.name.startswith("SHA256SUMS")}
(out / f"SHA256SUMS_{version}.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps({"source_files": len(paths), "artifacts": manifest}, indent=2))
