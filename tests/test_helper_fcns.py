from src.erik import load_mzml_data
from src.erik_constants import THOMAS_SAMPLES


class TestLoadMzmlData:
    @staticmethod
    def test_smoke():
        load_mzml_data(samples=THOMAS_SAMPLES[:2])
