"""Build a portable, checksummed snapshot of the completed OI-CBF work."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-11-08"
DEST = ROOT / "archive"
SUMMARY = DEST / f"{DATE}_工作归档.md"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def build(paper):
    DEST.mkdir(exist_ok=True)
    files = {}

    def add(path, name=None):
        files[name or path.relative_to(ROOT).as_posix()] = path.read_bytes()

    for name in ["README.md", "requirements.txt", ".gitignore", "hero.png"]:
        path = ROOT / name
        if path.is_file():
            add(path)
    for folder in ["double_integrator", "aircraft"]:
        for pattern in ["*.py", "*.m"]:
            for path in (ROOT / folder).rglob(pattern):
                add(path)
    for path in (ROOT / "matlab").rglob("*"):
        if path.is_file() and "results" not in path.relative_to(ROOT / "matlab").parts and path.suffix in [".m", ".md"]:
            add(path)
    for path in (ROOT / "tests").glob("*.py"):
        add(path)
    for path in (ROOT / "test_results").glob("*"):
        if path.is_file() and path.suffix in [".json", ".mat", ".npz", ".md"]:
            add(path)
    stale = {"validation_summary.mat", "code_analyzer.mat"}
    for path in (ROOT / "test_results" / "matlab").rglob("*"):
        relative = path.relative_to(ROOT / "test_results" / "matlab")
        if not path.is_file() or path.name in stale or any("entrypoint" in part for part in relative.parts):
            continue
        if path.suffix in [".json", ".mat", ".png", ".m"]:
            add(path)
    summary_source = SUMMARY if SUMMARY.is_file() else ROOT / "ARCHIVE_SUMMARY.md"
    add(summary_source, "ARCHIVE_SUMMARY.md")
    add(Path(__file__))
    add(paper, "paper/1048550arXiv251005436.pdf")

    # Verify original source remains equal to the hashes captured during audit.
    unchanged = True
    for model in ["double_integrator", "aircraft"]:
        audit = json.loads((ROOT / "test_results" / f"{model}_paper_audit.json").read_text(encoding="utf-8"))
        for path, expected in audit["source_sha256"].items():
            if digest((ROOT / path).read_bytes()) != expected:
                unchanged = False
    if not unchanged:
        raise RuntimeError("Original source differs from the audit snapshot; review before archiving.")

    result_dir = ROOT / "test_results" / "matlab"
    results = {name: json.loads((result_dir / f"{name}.json").read_text(encoding="utf-8")) for name in [
        "math_validation", "simulation_validation", "projection_validation", "paper_profile_validation",
        "toolbox_validation", "symbolic_validation", "parallel_validation", "code_analyzer"]}
    manifest = {
        "archive_date": DATE,
        "project": "OI-CBF",
        "source_unchanged_since_audit": unchanged,
        "file_count_excluding_manifest": len(files),
        "uncompressed_bytes_excluding_manifest": sum(map(len, files.values())),
        "validation_snapshot": results,
        "excluded": [".git", ".test-environment", "MATLAB installations", "stale aggregate MAT reports", "duplicate entrypoint output folders"],
        "files": [{"path": name, "bytes": len(data), "sha256": digest(data)} for name, data in sorted(files.items())],
    }
    manifest_data = (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    manifest_path = DEST / f"OI-CBF_{DATE}_归档清单.json"
    manifest_path.write_bytes(manifest_data)
    files["ARCHIVE_MANIFEST.json"] = manifest_data
    target = DEST / f"OI-CBF_{DATE}_完整归档.zip"
    prefix = f"OI-CBF_archive_{DATE}/"
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, data in sorted(files.items()):
            archive.writestr(prefix + name, data)
    with zipfile.ZipFile(target) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("ZIP integrity check failed.")
        for entry in manifest["files"]:
            if digest(archive.read(prefix + entry["path"])) != entry["sha256"]:
                raise RuntimeError(f"Manifest mismatch: {entry['path']}")
    checksum = digest(target.read_bytes())
    target.with_suffix(".zip.sha256").write_text(checksum + "  " + target.name + "\n", encoding="utf-8")
    print(json.dumps({"archive": str(target), "summary": str(SUMMARY), "files_including_manifest": len(files),
                      "archive_bytes": target.stat().st_size, "sha256": checksum,
                      "integrity": "ZIP CRC and every manifest SHA-256 verified"}, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    local_paper = ROOT / "paper" / "1048550arXiv251005436.pdf"
    default_paper = local_paper if local_paper.is_file() else Path(r"C:\Users\lu_ke\Downloads\1048550arXiv251005436.pdf")
    parser.add_argument("--paper", type=Path, default=default_paper)
    build(parser.parse_args().paper)
