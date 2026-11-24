import gzip
import json
import os
import tempfile
import unittest
from collections import Counter

import sf7_link_analysis as s

ADV_PK = "AB" * 32


def obs(observer, snr, path, radio="869.618,62.5,8,8"):
    return {"observer_id": observer, "snr": snr, "path_json": json.dumps(path), "radio": radio}


def rec(observations, payload_type=5, pubkey=None, route_type=1):
    r = {"payload_type": payload_type, "route_type": route_type, "observations": observations}
    if pubkey:
        key = "pubKey" if payload_type == 4 else "srcPubKey"
        r["decoded_json"] = json.dumps({key: pubkey.lower()})
    return r


class ExtractTests(unittest.TestCase):
    def test_last_hop_is_transmitter_and_bins(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", -7.6, ["1863", "84ca"])]))
        self.assertEqual(a.links[("84CA", "O1")], Counter({-7.5: 1}))
        self.assertEqual(a.adj[("1863", "84CA")], 1)

    def test_other_radio_skipped(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 1, ["AA"], radio="869.618,62.5,7,5")]))
        self.assertEqual(a.meta["radio_other"], 1)
        self.assertEqual(len(a.links), 0)

    def test_zero_hop_advert_uses_originator(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 3.0, [])], payload_type=4, pubkey=ADV_PK))
        self.assertIn((ADV_PK, "O1"), a.links)

    def test_radio_profile(self):
        self.assertEqual(s.radio_profile("869.617981,62.5,8,5"), "ok")
        self.assertEqual(s.radio_profile("869.6179809,62.5,8,8"), "ok")
        self.assertEqual(s.radio_profile("869.625,62.5,8,8"), "other")
        self.assertEqual(s.radio_profile("869.600,62.5,8,8"), "other")
        self.assertEqual(s.radio_profile("LoRa 869.618MHz"), "unknown")

    def test_direct_and_trace_paths_not_attributed_but_counted_as_receipt(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 5.0, ["AA", "BB"])], route_type=2))
        a.add_record(rec([obs("o2", 5.0, ["AA"])], payload_type=9, route_type=1))
        self.assertEqual(len(a.links), 0)
        self.assertEqual(a.meta["skip_direct_or_trace"], 2)
        self.assertEqual(a.meta["deliveries"], 2)

    def test_anon_req_src_pubkey(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 1.0, [])], payload_type=7, pubkey=ADV_PK, route_type=1))
        self.assertIn((ADV_PK, "O1"), a.links)
        b = s.Acc()  # empty DIRECT route = last forwarded copy, sender unknown
        b.add_record(rec([obs("o1", 1.0, [])], payload_type=7, pubkey=ADV_PK, route_type=2))
        self.assertEqual(len(b.links), 0)

    def test_tx_direction_and_bad_hops(self):
        a = s.Acc()
        o = obs("o1", 1.0, ["AA"])
        o["direction"] = "tx"
        a.add_record(rec([o, obs("o2", 1.0, ["ZZ"]), obs("o3", 1.0, ["AA", "BBCC"])]))
        self.assertEqual((a.meta["not_rx"], a.meta["bad_hop"], a.meta["mixed_hop_width"]), (1, 1, 1))
        self.assertEqual(len(a.links), 0)

    def test_receipt_counted_even_if_transmitter_unknown(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 9.0, []), obs("o1", -5.0, ["AA"])], payload_type=2))
        self.assertEqual(a.best_obs["O1"], Counter({9.0: 1}))
        self.assertEqual(a.links[("AA", "O1")], Counter({-5.0: 1}))

    def test_duplicate_observation_ids(self):
        a = s.Acc()
        seen = set()
        o = obs("o1", 1.0, ["AA"])
        o["id"] = 7
        a.add_record(rec([o]), seen)
        a.add_record(rec([dict(o)]), seen)
        self.assertEqual(a.meta["dup_observation"], 1)
        self.assertEqual(sum(a.links[("AA", "O1")].values()), 1)

    def test_zero_hop_unknown_originator_skipped(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", 3.0, [])], payload_type=2))
        self.assertEqual(a.meta["skip_unknown_tx"], 1)

    def test_best_copy_per_observer(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", -12, ["AA"]), obs("o1", -3, ["BB"]), obs("o2", -9, ["AA"])]))
        self.assertEqual(a.best, Counter({-3.0: 1, -9.0: 1}))
        self.assertEqual(a.meta["deliveries"], 2)
        self.assertEqual(a.best_obs["O1"], Counter({-3.0: 1}))

    def test_json_roundtrip_and_file(self):
        a = s.Acc()
        a.add_record(rec([obs("o1", -12, ["AA", "BB"])]))
        b = s.Acc.from_json(json.loads(json.dumps(a.to_json())))
        self.assertEqual(b.links, a.links)
        self.assertEqual(b.adj, a.adj)
        self.assertEqual(b.best_obs, a.best_obs)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.jsonl.gz")
            with gzip.open(p, "wt") as fh:
                fh.write(json.dumps(rec([obs("o1", -1, ["CC"])])) + "\n\nnot json\n")
            c = s.extract_file(p)
            self.assertEqual(c.meta["used"], 1)
            self.assertEqual(c.meta["bad_line"], 1)


class StatsTests(unittest.TestCase):
    def test_quantile_and_frac(self):
        h = Counter({-10.0: 1, -5.0: 1, 0.0: 2})
        self.assertEqual(s.hist_quantile(h, 0.25), -10.0)
        self.assertEqual(s.hist_quantile(h, 0.5), -5.0)
        self.assertAlmostEqual(s.frac_below(h, -5.0), 0.25)
        self.assertIsNone(s.frac_below(Counter(), 0))

    def test_resolve_prefix(self):
        idx = {"1863AA": {"lat": 50.0, "lon": 7.0, "role": "repeater"},
               "84CA00": {"lat": 50.0, "lon": 7.0, "role": "repeater"},
               "84CB00": {"lat": 52.0, "lon": 9.0, "role": "repeater"},
               "OBS": {"lat": 50.01, "lon": 7.0}}
        keys = sorted(idx)
        self.assertEqual(s.resolve_prefix("1863", "OBS", idx, keys, {}), "1863AA")
        self.assertEqual(s.resolve_prefix("84C", "OBS", idx, keys, {}), "84CA00")
        idx["84CC00"] = {"lat": 50.0, "lon": 7.01, "role": "repeater"}
        keys = sorted(idx)
        self.assertIsNone(s.resolve_prefix("84C", "OBS", idx, keys, {}))
        # a candidate without position keeps the prefix ambiguous
        del idx["84CC00"]
        idx["84CD00"] = {"lat": None, "lon": None, "role": "repeater"}
        keys = sorted(idx)
        self.assertIsNone(s.resolve_prefix("84C", "OBS", idx, keys, {}))
        self.assertIsNone(s.resolve_prefix("FFFF", "OBS", idx, keys, {}))
        self.assertIsNone(s.resolve_prefix("84C", "OBS", {"84CA00": idx["84CA00"], "84CB00": idx["84CB00"],
                                                          "OBS": idx["OBS"]}, ["84CA00", "84CB00", "OBS"], {}, strict=True))

    def test_analyse_floor_and_loss(self):
        a = s.Acc()
        # 100 obs: 10 at -12 (below floor+2.5 when floor is -12), 90 at +5
        a.links[("AA", "O1")] = Counter({-12.0: 10, 5.0: 90})
        a.best = Counter({-12.0: 10, 5.0: 90})
        a.best_obs["O1"] = Counter({-12.0: 10, 5.0: 90})
        idx = {"AA11": {"lat": 50.0, "lon": 7.0, "role": "repeater"},
               "O1": {"lat": 50.0, "lon": 7.05, "observer": True}}
        out = s.analyse(a, idx, center=(50.0, 7.0), radii=[10], shifts=[2.5], floor_q=0.05)
        self.assertEqual(out["sf8_floor"], -12.0)
        res = out["shifts"]["2.5"]
        self.assertEqual(res["sf7_floor"], -9.5)
        self.assertAlmostEqual(res["delivery_loss"], 0.10)
        area = res["areas"]["10 km"]
        self.assertEqual(area["observations"], 100)
        self.assertAlmostEqual(area["delivery_loss"], 0.10)
        self.assertEqual(area["link_classes"], {"degraded": 1})
        self.assertEqual(area["observer_best_link"], {"ok": 1})


class NeighbourTests(unittest.TestCase):
    def test_dedup_same_advert_across_polls(self):
        lines = [
            json.dumps({"ts": "2026-11-23T10:00:00+00:00", "type": "neighbours",
                        "neighbours": [{"id": "d0beefa6", "secs_ago": 100, "snr_db": 10.75}]}),
            # same advert 10 min later (heard time identical; in another file -> global dedup)
            json.dumps({"ts": "2026-11-23T10:10:00+00:00", "type": "neighbours",
                        "neighbours": [{"id": "D0BEEFA6", "secs_ago": 700, "snr_db": 10.75}]}),
            # a new advert
            json.dumps({"ts": "2026-11-23T14:10:00+00:00", "type": "neighbours",
                        "neighbours": [{"id": "D0BEEFA6", "secs_ago": 5, "snr_db": 9.0}]}),
            json.dumps({"ts": "2026-11-23T14:10:00+00:00", "type": "stats"}),
        ]
        per = s.neighbour_samples(lines)
        self.assertEqual(per["D0BEEFA6"], [10.75, 9.0])

    def test_distinct_adverts_close_in_time_kept(self):
        lines = [json.dumps({"ts": "2026-11-23T10:00:00+00:00", "type": "neighbours",
                             "neighbours": [{"id": "AA", "secs_ago": 70, "snr_db": 1.0}]}),
                 json.dumps({"ts": "2026-11-23T10:01:00+00:00", "type": "neighbours",
                             "neighbours": [{"id": "AA", "secs_ago": 5, "snr_db": 2.0}]})]
        self.assertEqual(s.neighbour_samples(lines)["AA"], [1.0, 2.0])

    def test_since_filter(self):
        from datetime import datetime
        lines = [json.dumps({"ts": "2026-10-27T00:00:30+00:00", "type": "neighbours",
                             "neighbours": [{"id": "AA", "secs_ago": 60, "snr_db": 1}]})]  # heard before cutoff
        self.assertEqual(s.neighbour_samples(lines, datetime.fromisoformat("2026-10-27T00:00:00+00:00")), {})


if __name__ == "__main__":
    unittest.main()
