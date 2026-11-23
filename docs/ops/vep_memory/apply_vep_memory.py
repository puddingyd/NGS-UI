#!/usr/bin/env python3
"""Preview / install a lower VEP_ANNOTATE memory request in one profile.

Why: in the dgm profile ``process_high`` asks for 32 cpus and 64 GB. VEP is
a process_high task but really peaks at ~7 GB (9 DRAGEN WGS, 2026-11-22),
and on the 125 GB server two 64 GB requests do not fit, so Nextflow ran
every VEP one after another (~7.5 h for 9 WGS). The memory directive is only
a scheduling request for the local executor (Apptainer gets no memory
limit), and VEP's script uses ``task.cpus`` but not ``task.memory``, so
lowering it changes neither the results nor the Nextflow cache.

One line is added inside ``profiles { <profile> { process { ... } } }``
right after its ``withLabel: process_high`` line::

    withName: 'VEP_ANNOTATE' { memory = 16.GB }

(``withName`` takes priority over ``withLabel``; cpus stay at the label's
value.) Default mode only PREVIEWS: the edited copy is written to
--stage-dir and the diff is printed. ``--install`` checks the config did not
change since the preview, writes a timestamped backup next to it and prints
the restore command.
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

PROFILES = ("local", "dgm", "dgx")


def die(message):
    print("ERROR: " + message, file=sys.stderr)
    sys.exit(1)


def sha256(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def read(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return handle.read()


def edit_config(text, profile, memory):
    if re.search(r"withName\s*:\s*['\"]?VEP_ANNOTATE\b", text):
        die("the config already has a withName: 'VEP_ANNOTATE' selector (already installed?)")
    lines = text.splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if re.match(r"^\s*%s\s*\{\s*$" % re.escape(profile), line)]
    if len(starts) != 1:
        die("expected exactly one '%s {' profile block, found %d" % (profile, len(starts)))
    start = starts[0]
    others = re.compile(r"^\s*(%s)\s*\{\s*$" % "|".join(p for p in PROFILES if p != profile))
    target = None
    for i in range(start + 1, len(lines)):
        if others.match(lines[i]):
            break
        if re.match(r"^\s*withLabel\s*:\s*['\"]?process_high['\"]?\s*\{", lines[i]):
            target = i
            break
    if target is None:
        die("no 'withLabel: process_high' line inside the '%s' profile" % profile)
    indent = re.match(r"^(\s*)", lines[target]).group(1)
    newline = "\r\n" if lines[target].endswith("\r\n") else "\n"
    added = (
        "{i}// VEP 實際約 7 GB（9 例 WGS 實測）；只降排程宣告，讓兩個 VEP 可同時跑。"
        "不影響結果與 -resume cache（VEP 只用 task.cpus）{n}"
        "{i}withName: 'VEP_ANNOTATE' {{ memory = {m} }}{n}"
    ).format(i=indent, n=newline, m=memory)
    return "".join(lines[:target + 1]) + added + "".join(lines[target + 1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--config", default="/home/pipeline/tertiary_code/nextflow_tertiary.config")
    parser.add_argument("--profile", default="dgm", choices=PROFILES)
    parser.add_argument("--memory", default="16.GB")
    parser.add_argument("--stage-dir", default=os.path.expanduser("~/vep_memory_deploy"))
    parser.add_argument("--install", action="store_true")
    args = parser.parse_args()

    if not re.match(r"^\d+\.(GB|MB)$", args.memory):
        die("--memory must look like 16.GB")
    if not os.path.isfile(args.config):
        die("not found: " + args.config)
    os.makedirs(args.stage_dir, exist_ok=True)
    staged = os.path.join(args.stage_dir, os.path.basename(args.config))
    manifest = os.path.join(args.stage_dir, "manifest.json")

    if not args.install:
        new_text = edit_config(read(args.config), args.profile, args.memory)
        with open(staged, "w", encoding="utf-8", newline="") as handle:
            handle.write(new_text)
        with open(manifest, "w", encoding="utf-8") as handle:
            json.dump({"config": args.config, "sha256": sha256(args.config)}, handle)
        print("=" * 72)
        print("PREVIEW ONLY — the config was NOT changed. Staged copy: " + staged)
        print("=" * 72)
        sys.stdout.writelines(difflib.unified_diff(
            read(args.config).splitlines(keepends=True), new_text.splitlines(keepends=True),
            fromfile=args.config, tofile=args.config + " (new)", n=3))
        return 0

    if not (os.path.isfile(staged) and os.path.isfile(manifest)):
        die("nothing staged; run without --install first and review the diff")
    with open(manifest, encoding="utf-8") as handle:
        recorded = json.load(handle)
    if recorded.get("config") != args.config or recorded.get("sha256") != sha256(args.config):
        die("the config changed since the preview; run the preview again")
    backup = "{}.bak-{}".format(args.config, datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    shutil.copy2(args.config, backup)
    shutil.copyfile(staged, args.config)
    os.remove(manifest)
    print("installed {}  (backup: {})".format(args.config, backup))
    print("\nTo undo, run:\n  cp -p {} {}".format(backup, args.config))
    return 0


if __name__ == "__main__":
    sys.exit(main())
