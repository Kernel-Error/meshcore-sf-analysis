import unittest

from lora_airtime import airtime_s, meshcore_preamble, sensitivity_dbm
import lora_airtime as la


class LoraAirtimeTest(unittest.TestCase):
    def test_reference_sf7_125k(self):
        # Semtech LoRa calculator: SF7/125 kHz/CR4-5, 8 preamble, 10 B, CRC on -> 41.2 ms
        self.assertAlmostEqual(airtime_s(10, 7, 125000, 1, preamble=8) * 1000, 41.2, places=1)

    def test_low_data_rate_optimisation_kicks_in(self):
        # SF12/125 kHz symbol time 32.8 ms -> DE=1, longer than without
        self.assertGreater(airtime_s(10, 12, 125000, 1), 1.0)

    def test_higher_coding_rate_costs_airtime(self):
        self.assertGreater(airtime_s(50, 8, 62500, 4), airtime_s(50, 8, 62500, 1))

    def test_meshcore_preamble_rule(self):
        self.assertEqual(meshcore_preamble(8), 32)
        self.assertEqual(meshcore_preamble(9), 16)

    def test_default_uses_meshcore_preamble(self):
        # SF8/62.5k: 16 extra preamble symbols = 16 * 4.096 ms = ~65.5 ms
        diff = airtime_s(20, 8, 62500, 4) - airtime_s(20, 8, 62500, 4, preamble=16)
        self.assertAlmostEqual(diff * 1000, 65.5, places=1)

    def test_published_table_values(self):
        # Part 1 table (corrected 2026-11-04): MeshCore SF8/62.5k/CR4-8 and LongFast SF11/250k/CR4-5 (preamble 16)
        self.assertAlmostEqual(airtime_s(20, 8, 62500, 4) * 1000, 378, delta=1)
        self.assertAlmostEqual(airtime_s(180, 8, 62500, 4) * 1000, 1689, delta=1)
        self.assertAlmostEqual(airtime_s(180, 11, 250000, 1, preamble=16) * 1000, 1583, delta=1)

    def test_sensitivity(self):
        self.assertAlmostEqual(sensitivity_dbm(8, 62500), -130.0, places=1)


class Sx126xDriverParityTests(unittest.TestCase):
    """Values cross-checked against Semtech sx126x_get_lora_time_on_air_numerator (integer ms ceil)."""

    def _driver_ms(self, pl, sf, bw, cr, pre, crc=True, implicit=False, ldro=False):
        num = (pl << 3) + (16 if crc else 0) - 4 * sf + (0 if implicit else 20)
        if sf <= 6:
            den = 4 * sf
        else:
            num += 8
            den = 4 * (sf - 2) if ldro else 4 * sf
        num = max(num, 0)
        inter = ((num + den - 1) // den) * (cr + 4) + pre + 12
        if sf <= 6:
            inter += 2
        n = (4 * inter + 1) * (1 << (sf - 2))
        return -(-1000 * n // bw)

    def test_parity(self):
        for pl in (10, 40, 100, 200):
            for sf in (5, 6, 7, 8, 9):
                for cr in (1, 4):
                    pre = la.meshcore_preamble(sf)
                    got = la.airtime_s(pl, sf, 62500, cr) * 1000
                    self.assertAlmostEqual(got, self._driver_ms(pl, sf, 62500, cr, pre), delta=1.0,
                                           msg=f"pl={pl} sf={sf} cr={cr}")


if __name__ == "__main__":
    unittest.main()
