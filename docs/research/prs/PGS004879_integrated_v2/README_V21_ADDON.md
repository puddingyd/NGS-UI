# PGS004879 integrated v2.1 — additive execution layer

This is an **add-on**, not a replacement for `PGS004879_integrated_v2`.

It installs new `pipeline_v21.py`, `v21_raw_resume.py`, `tests/test_v21.py`,
`V21_SOP.md`, manifests, and a 10-file baseline SHA256 manifest into the
existing v2 folder. It refuses to overwrite files, and verifies the v2 source
components before installing. No changes to genotype or PRS scoring code.

**READ:** `V21_SOP.md` (Traditional Chinese).

**Recommended acceptance:** (1) unit tests and model lock; (2) 3-case reuse
regression; (3) optional 3-case fresh gVCF regression; (4) separate failure
injection and concurrency tests before scale-up.

**Current validation:** 25 synthetic/mock unit tests PASS in the artifact
creation environment. No v2.1 job has been executed on the user's server.
Original v2 three-case server validation remains independently established.
