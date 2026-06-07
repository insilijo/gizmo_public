"""Deploy a built GIZMO bundle to GitHub Release.

Works against private GitHub repos via `gh` CLI auth. Zenodo auto-mirror
is paused while the repo is private — once made public, enable the
webhook at https://zenodo.org/account/settings/github and the next
release will auto-DOI.

Workflow:
  1. Tar the bundle directory into bundle-v{X.Y.Z}.tar.gz
  2. Compute sha256 + size of the tarball
  3. Create a git tag for the version (e.g., bundle-v1.0.0)
  4. Push the tag to origin
  5. Create a GitHub Release attached to the tag via `gh release create`
  6. Upload the tarball as a release asset

Usage::

    # Tag and push the current human_full bundle as v1.0.0:
    python -m gizmo.build.deploy_bundle data/processed/human_full

    # Use a custom tag name (default: bundle-v{version} from manifest):
    python -m gizmo.build.deploy_bundle data/processed/human_full --tag bundle-v1.0.0

    # Dry run — build the tar but don't push or release:
    python -m gizmo.build.deploy_bundle data/processed/human_full --dry-run

Prerequisites:
  - `gh` (GitHub CLI) installed + authenticated (`gh auth login`)
  - Repo has a configured remote (origin)

When ready for paper-grade DOI:
  - Either make the repo public + enable Zenodo webhook (auto-mirrors)
  - Or manually upload the tarball to zenodo.org for a direct DOI
"""
from __future__ import annotations

import argparse, json, subprocess, sys, hashlib, tarfile, shutil
from pathlib import Path
from datetime import datetime, timezone


def _sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf: break
            h.update(buf)
    return h.hexdigest()


def package_bundle(bundle_dir: Path, out_dir: Path) -> Path:
    """Tar.gz the bundle directory; return path to the tarball."""
    manifest_path = bundle_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"manifest.json not found in {bundle_dir}")
    manifest = json.loads(manifest_path.read_text())
    version = manifest.get("bundle_version", "unversioned")
    name = manifest.get("graph_name", bundle_dir.name)

    out_dir.mkdir(parents=True, exist_ok=True)
    tarball = out_dir / f"{name}-v{version}.tar.gz"
    if tarball.exists(): tarball.unlink()

    print(f"  packaging {bundle_dir} → {tarball}")
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(bundle_dir, arcname=name)
    return tarball


def gh(cmd: list[str], dry_run: bool = False) -> str:
    full = ["gh"] + cmd
    print(f"  $ {' '.join(full)}")
    if dry_run: return "(dry run)"
    r = subprocess.run(full, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"gh command failed: {' '.join(full)}")
    return r.stdout.strip()


def git(cmd: list[str], dry_run: bool = False, cwd: Path | None = None) -> str:
    full = ["git"] + cmd
    print(f"  $ {' '.join(full)}")
    if dry_run: return "(dry run)"
    r = subprocess.run(full, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"git command failed: {' '.join(full)}")
    return r.stdout.strip()


def deploy(bundle_dir: Path, tag: str | None = None, dry_run: bool = False,
            release_notes: str | None = None) -> dict:
    """Tar, tag, push, create GitHub release with the tarball attached."""
    bundle_dir = Path(bundle_dir)
    if not bundle_dir.exists():
        raise FileNotFoundError(bundle_dir)

    manifest = json.loads((bundle_dir / "manifest.json").read_text())
    version = manifest.get("bundle_version", "unversioned")
    name = manifest.get("graph_name", bundle_dir.name)

    # Default tag: bundle-v{version}
    if tag is None:
        tag = f"bundle-v{version}"

    print(f"=== Deploying bundle {name} v{version} as tag {tag} ===")

    # 1. Package
    out_dir = bundle_dir.parent / "_release_artifacts"
    tarball = package_bundle(bundle_dir, out_dir)
    sha = _sha256_file(tarball)
    size = tarball.stat().st_size
    print(f"  tarball: {size/1024/1024:.1f}MB  sha256: {sha[:16]}…")

    # Default release notes
    if release_notes is None:
        nodes = (manifest.get("node_counts") or {})
        sources = [s.get("name") if isinstance(s, dict) else s
                   for s in (manifest.get("sources") or [])]
        release_notes = (
            f"# {name} v{version}\n\n"
            f"Multi-source GIZMO graph bundle.\n\n"
            f"## Build metadata\n"
            f"- Built at: {manifest.get('built_at')}\n"
            f"- Versioned at: {manifest.get('versioned_at')}\n"
            f"- Sources: {', '.join(sources)}\n\n"
            f"## Artifact\n"
            f"- File: `{tarball.name}`\n"
            f"- Size: {size:,} bytes\n"
            f"- SHA-256: `{sha}`\n\n"
            f"## Verify after download\n"
            f"```bash\n"
            f"tar -xzf {tarball.name}\n"
            f"python -m gizmo.build.build_human_full verify {name}\n"
            f"```\n"
        )

    # 2. Git tag + push
    git(["tag", "-a", tag, "-m", f"GIZMO bundle {name} v{version}"], dry_run=dry_run)
    git(["push", "origin", tag], dry_run=dry_run)

    # 3. GitHub release
    notes_file = out_dir / "_release_notes.md"
    notes_file.write_text(release_notes)
    gh(["release", "create", tag,
        str(tarball),
        "--title", f"{name} v{version}",
        "--notes-file", str(notes_file)],
       dry_run=dry_run)

    # 4. Get release URL
    if not dry_run:
        repo_url = git(["remote", "get-url", "origin"], dry_run=False)
        repo_url = (repo_url
                     .replace("git@github.com:", "https://github.com/")
                     .replace(".git", ""))
        release_url = f"{repo_url}/releases/tag/{tag}"
    else:
        release_url = "(dry run)"

    return {
        "tag": tag, "version": version, "name": name,
        "tarball": str(tarball),
        "tarball_size_bytes": size,
        "tarball_sha256": sha,
        "release_url": release_url,
        "zenodo_status": (
            "Repo is private; Zenodo auto-mirror is paused. When ready, "
            "either (a) make the repo public + enable the Zenodo-GitHub "
            "webhook and re-deploy, or (b) manually upload the tarball "
            "to zenodo.org for a direct DOI."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description="Deploy a GIZMO bundle to GitHub Release")
    parser.add_argument("bundle_dir", help="Path to bundle directory")
    parser.add_argument("--tag", default=None, help="Override default tag (bundle-v{version})")
    parser.add_argument("--dry-run", action="store_true", help="Don't push or upload")
    parser.add_argument("--notes-file", default=None,
                          help="Path to custom release notes markdown")

    args = parser.parse_args()
    notes = None
    if args.notes_file:
        notes = Path(args.notes_file).read_text()
    result = deploy(Path(args.bundle_dir), tag=args.tag,
                     dry_run=args.dry_run, release_notes=notes)
    print()
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
