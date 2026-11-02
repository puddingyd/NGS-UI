import argparse
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pipeline_v21 as v21


class BatchSuccessRegression(unittest.TestCase):

    def test_three_namespaces_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "sample\tmode\tgvcf\traw_source\n"
                + "".join(
                    f"26T0097{i}\treuse\t\t/raw/{i}\n"
                    for i in range(3)
                )
            )

            args = argparse.Namespace(
                config="/unused",
                manifest=str(manifest),
                run_id="testbatch",
                jobs=1,
                resume=False,
                with_legacy_ancestry=False,
            )

            records = [
                {
                    "sample": f"26T0097{i}",
                    "execution_status": "PASS",
                    "model": model,
                }
                for i in range(3)
                for model in ("A", "B")
            ]

            with (
                patch.object(
                    v21.v1, "cfg_read",
                    return_value={"project": tmp}
                ),
                patch.object(v21, "run_one", return_value=None)
                as runner,
                patch.object(
                    v21, "summarize", return_value=records
                ),
            ):
                v21.batch(args)

            self.assertEqual(runner.call_count, 3)

    def test_failure_reports_unfinished(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            manifest.write_text(
                "sample\tmode\tgvcf\traw_source\n"
                "26T00970\treuse\t\t/raw/00970\n"
            )

            args = argparse.Namespace(
                config="/unused",
                manifest=str(manifest),
                run_id="testbatch",
                jobs=1,
                resume=False,
                with_legacy_ancestry=False,
            )

            with (
                patch.object(
                    v21.v1, "cfg_read",
                    return_value={"project": tmp}
                ),
                patch.object(
                    v21, "run_one",
                    side_effect=RuntimeError("synthetic failure")
                ),
                patch.object(
                    v21, "summarize",
                    return_value=[]
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "Batch unfinished"
                ):
                    v21.batch(args)


if __name__ == "__main__":
    unittest.main()
