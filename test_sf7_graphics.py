import json
import math
import unittest
from collections import Counter

import sf7_graphics as g
import sf7_link_analysis as sla


class ConfigTests(unittest.TestCase):
    def test_reference_and_ordering(self):
        rows = {r["label"]: r for r in g.config_table()}
        self.assertEqual(rows["SF8/CR8"]["airtime_rel"], 1.0)
        self.assertEqual(rows["SF8/CR8"]["budget_vs_ref_db"], 0.0)
        self.assertAlmostEqual(rows["SF7/CR8"]["budget_vs_ref_db"], -2.5, places=1)
        self.assertAlmostEqual(rows["SF9/CR8"]["budget_vs_ref_db"], 2.5, places=1)
        self.assertLess(rows["SF7/CR5"]["airtime_ms"], rows["SF8/CR5"]["airtime_ms"])
        self.assertGreater(rows["SF9/CR8"]["airtime_rel"], 1.5)
        self.assertAlmostEqual(rows["SF8/CR5 125k"]["budget_vs_ref_db"], -3.0, places=1)
        self.assertEqual(rows["SF8/CR5 125k"]["note"], "band")
        self.assertEqual(rows["SF7/CR5 31k"]["note"], "clock")
        self.assertAlmostEqual(rows["SF6/CR5"]["airtime_ms"], 221.4, delta=0.5)

    def test_aloha(self):
        self.assertEqual(g.aloha_success(0), 1.0)
        self.assertAlmostEqual(g.aloha_success(0.5), math.exp(-1))

    def test_de_num(self):
        self.assertEqual(g.de_num(-7.5), "−7,5")
        self.assertEqual(g.de_num(12.0, 0), "12")


class LoadTests(unittest.TestCase):
    def _st(self, ts, up, rx, tx):
        return json.dumps({"ts": ts, "type": "stats", "core": {"uptime_secs": up},
                           "radio": {"rx_air_secs": rx, "tx_air_secs": tx}})

    def _day(self, date="2026-11-20"):
        return [self._st(f"{date}T{i // 60:02d}:{i % 60:02d}:00+00:00", 1000 + 60 * i, 18 * i, 3 * i)
                for i in range(0, 24 * 60)]

    def test_shares_sorted_dedup_and_reboot(self):
        lines = self._day()
        lines = lines[::-1] + lines[:5]  # reversed order + duplicates
        lines.append(self._st("2026-11-20T23:59:30+00:00", 10, 0, 0))  # reboot -> no delta
        out = g.channel_load(lines, min_hours=23.0)
        d = out["2026-11-20"]
        self.assertAlmostEqual(d["rx_share"], 0.3)
        self.assertAlmostEqual(d["tx_share"], 0.05)

    def test_partial_day_dropped_and_window(self):
        from datetime import datetime
        lines = self._day()
        self.assertEqual(g.channel_load(lines, until=datetime.fromisoformat("2026-11-20T12:00:00+00:00")), {})

    def test_wall_time_mismatch_skipped(self):
        lines = [self._st("2026-11-20T00:00:00+00:00", 0, 0, 0), self._st("2026-11-20T00:01:00+00:00", 600, 60, 0)]
        self.assertEqual(g.channel_load(lines, min_hours=0), {})


class AreaTests(unittest.TestCase):
    def test_area_filter_and_equal_observer(self):
        a = sla.Acc()
        a.links[("AA", "NEAR")] = Counter({-12.0: 1, 0.0: 3})
        a.links[("AA", "FAR")] = Counter({-12.0: 4})
        a.best_obs["NEAR"] = Counter({0.0: 4})
        a.best_obs["FAR"] = Counter({-12.0: 4})
        idx = {"NEAR": {"lat": 50.58, "lon": 6.95}, "FAR": {"lat": 52.5, "lon": 13.4}}
        oh, dh, used = g.area_hists(a, idx, g.CENTER, 15)
        self.assertEqual(sum(oh.values()), 4)
        self.assertEqual(dh, Counter({0.0: 4}))
        self.assertEqual(used, {"NEAR"})
        oh, dh, used = g.area_hists(a, idx, g.CENTER, None)
        self.assertEqual(sum(oh.values()), 8)
        self.assertEqual(g.loss_curve(dh, [-11.0])[0][1], 0.5)
        eq = g.equal_observer_exposure(a, used, -11.0)
        self.assertEqual((eq["mean"], eq["min"], eq["max"], eq["n"]), (0.5, 0.0, 1.0, 2))
        self.assertIsNone(g.equal_observer_exposure(a, set(), -11.0))

    def test_cdf_strictly_below(self):
        h = Counter({-10.0: 1, -5.0: 1})
        self.assertEqual(g.cdf_below(h, [-10.0, -9.75, -5.0, -4.75]), [0.0, 50.0, 50.0, 100.0])
        self.assertEqual(g.cdf_below(Counter(), [0.0]), [])


class DistanceTests(unittest.TestCase):
    def test_classes_and_share(self):
        rows = [{"km": 1, "n": 100, "loss": 0.0}, {"km": 5, "n": 100, "loss": 0.2},
                {"km": 5, "n": 300, "loss": 0.5}, {"km": 5, "n": 5, "loss": 1.0}, {"km": 90, "n": 50, "loss": 1.0}]
        out = g.distance_classes(rows)
        self.assertEqual(out[0]["ok"], 1)
        self.assertEqual((out[1]["links"], out[1]["degraded"], out[1]["lost"]), (2, 1, 1))  # 0.5 -> '>= 50 %'
        self.assertAlmostEqual(out[1]["obs_share"], (20 + 150) / 400)
        self.assertIsNone(out[4]["obs_share"])
        self.assertEqual(g.share_ge10(rows), {"links": 2, "ge10": 2, "share": 1.0})


if __name__ == "__main__":
    unittest.main()
