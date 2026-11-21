#!/usr/bin/env python3
"""Stage, review and install the Pangolin score cache in the production pipeline.

Two edits are made to the production tertiary pipeline:

  modules/snv_annotation.nf   process PANGOLIN_SCORE replaced by PANGOLIN_SCORE.nf
  nextflow_tertiary.config    one line added after `pangolin_db = ...`:
                              pangolin_cache_db = "<cache path>"

Default mode only PREVIEWS: the edited copies are written to --stage-dir and a
diff against the production files is printed. Production files are untouched.

  --install   copies the staged files into place, after
              * checking the production files did not change since the preview
              * saving timestamped backups next to them (*.bak-YYYYmmdd-HHMMSS)
              It prints the exact commands to undo the installation.

Refuses to do anything when the files do not look as expected (process not
found exactly once, cache already installed, ...).
"""
import argparse
import datetime
import difflib
import hashlib
import json
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NEW_PROCESS_FILE = os.path.join(HERE, "PANGOLIN_SCORE.nf")
DEFAULT_CACHE = "${params.ref_dir}/tertiary/pangolin/cache/pangolin_scores.sqlite"


def die(message):
    print("ERROR: " + message, file=sys.stderr)
    sys.exit(1)


def sha256(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def read(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return handle.read()


def new_process_text():
    text = read(NEW_PROCESS_FILE)
    start = text.find("process PANGOLIN_SCORE {")
    if start < 0:
        die("PANGOLIN_SCORE.nf does not contain the process")
    return text[start:].rstrip("\n") + "\n"


def edit_module(text):
    lines = text.splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if line.startswith("process PANGOLIN_SCORE {")]
    if len(starts) != 1:
        die("expected exactly one 'process PANGOLIN_SCORE {' line in snv_annotation.nf, "
            "found {}".format(len(starts)))
    start = starts[0]
    end = next((i for i in range(start + 1, len(lines)) if lines[i].rstrip("\r\n") == "}"), None)
    if end is None:
        die("could not find the closing '}' of process PANGOLIN_SCORE")
    old_block = "".join(lines[start:end + 1])
    if "pangolin_cache.py" in old_block:
        die("PANGOLIN_SCORE already uses pangolin_cache.py (already installed?)")
    for marker in ("pangolin \\\\", "splice_for_pangolin.vcf", "02_pangolin"):
        if marker not in old_block:
            die("the current PANGOLIN_SCORE does not look like the expected v3.6 process "
                "(missing {!r}); not touching it".format(marker))
    if re.search(r"^\s*(?:process|workflow)\s", "".join(lines[start + 1:end]), re.M):
        die("unexpected nested process/workflow inside PANGOLIN_SCORE block; not touching it")
    return "".join(lines[:start]) + new_process_text() + "".join(lines[end + 1:])


def edit_config(text, cache_expr):
    if re.search(r"^\s*pangolin_cache_db\s*=", text, re.M):
        die("nextflow_tertiary.config already declares pangolin_cache_db")
    lines = text.splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if re.match(r"^\s*pangolin_db\s*=", line)]
    if len(hits) != 1:
        die("expected exactly one 'pangolin_db = ...' line in nextflow_tertiary.config, "
            "found {}".format(len(hits)))
    i = hits[0]
    indent = re.match(r"^(\s*)", lines[i]).group(1)
    newline = "\r\n" if lines[i].endswith("\r\n") else "\n"
    added = '{}pangolin_cache_db = "{}"{}'.format(indent, cache_expr, newline)
    comment = ("{}// Pangolin 分數快取（docs/ops/pangolin_cache）；設成 '' 即停用，"
               "行為與舊版相同{}").format(indent, newline)
    return "".join(lines[:i + 1]) + comment + added + "".join(lines[i + 1:])


def show_diff(old_path, new_text):
    diff = difflib.unified_diff(
        read(old_path).splitlines(keepends=True), new_text.splitlines(keepends=True),
        fromfile=old_path, tofile=old_path + " (new)", n=2)
    sys.stdout.writelines(diff)
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tertiary-code", default="/home/pipeline/tertiary_code")
    parser.add_argument("--stage-dir", default=os.path.expanduser("~/pangolin_cache_deploy"))
    parser.add_argument("--cache-expr", default=DEFAULT_CACHE,
                        help="value for pangolin_cache_db (default: %(default)s)")
    parser.add_argument("--install", action="store_true")
    args = parser.parse_args()

    module = os.path.join(args.tertiary_code, "modules", "snv_annotation.nf")
    config = os.path.join(args.tertiary_code, "nextflow_tertiary.config")
    script = os.path.join(args.tertiary_code, "scripts", "pangolin_cache.py")
    for path in (module, config):
        if not os.path.isfile(path):
            die("not found: " + path)
    os.makedirs(args.stage_dir, exist_ok=True)
    staged_module = os.path.join(args.stage_dir, "snv_annotation.nf")
    staged_config = os.path.join(args.stage_dir, "nextflow_tertiary.config")
    manifest_path = os.path.join(args.stage_dir, "manifest.json")

    if not args.install:
        new_module = edit_module(read(module))
        new_config = edit_config(read(config), args.cache_expr)
        with open(staged_module, "w", encoding="utf-8", newline="") as handle:
            handle.write(new_module)
        with open(staged_config, "w", encoding="utf-8", newline="") as handle:
            handle.write(new_config)
        with open(manifest_path, "w", encoding="utf-8") as handle:
            json.dump({module: sha256(module), config: sha256(config)}, handle, indent=1)
        print("=" * 72)
        print("PREVIEW ONLY — production files were NOT changed.")
        print("Staged copies: {}".format(args.stage_dir))
        print("=" * 72)
        show_diff(config, new_config)
        show_diff(module, new_module)
        if not os.path.isfile(script):
            print("WARNING: {} not found — copy pangolin_cache.py there before "
                  "installing.".format(script))
        else:
            ours = os.path.join(HERE, "pangolin_cache.py")
            same = sha256(script) == sha256(ours)
            print("pangolin_cache.py in scripts/: {}".format(
                "identical to the repo copy" if same else "DIFFERS from the repo copy — copy it again"))
        return 0

    if not (os.path.isfile(staged_module) and os.path.isfile(staged_config)
            and os.path.isfile(manifest_path)):
        die("nothing staged; run without --install first and review the diff")
    with open(manifest_path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    for path in (module, config):
        if manifest.get(path) != sha256(path):
            die("{} changed since the preview; run the preview again".format(path))
    if not os.path.isfile(script):
        die("{} not found; copy pangolin_cache.py there first".format(script))
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    restore = []
    for path, staged in ((module, staged_module), (config, staged_config)):
        backup = "{}.bak-{}".format(path, stamp)
        shutil.copy2(path, backup)
        shutil.copyfile(staged, path)
        restore.append("cp -p {} {}".format(backup, path))
        print("installed {}  (backup: {})".format(path, backup))
    os.remove(manifest_path)
    print("\nTo undo, run:")
    for line in restore:
        print("  " + line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
