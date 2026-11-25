"""Run the production script's embedded offline reference tests."""
import json
import subprocess
import sys
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "mapbiomas_viirs_collection1_final_1to1.py"


class ProductionSelfTest(unittest.TestCase):
    def test_reference_and_configuration(self):
        result = subprocess.run([sys.executable, str(SCRIPT), "--self-test"],
                                capture_output=True, text=True, check=True)
        report, _ = json.JSONDecoder().raw_decode(result.stdout.lstrip())
        self.assertTrue(report["passed"])
        self.assertTrue(all(report["configuration"].values()))
        self.assertTrue(report["numerical_reference"]["passed"])


if __name__ == "__main__":
    unittest.main()
