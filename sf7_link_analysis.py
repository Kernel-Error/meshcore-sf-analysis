#!/usr/bin/env python3
"""Estimate what an SF8 -> SF7 switch would mean for links seen by CoreScope.

Two steps:

  extract <dump.jsonl.gz>... --out links.json.gz [--workers N]
      Streams CoreScope 'packets-expanded' day dumps (stdlib only, runs on the
      Pi). Keeps only observations on 869.618 MHz / 62.5 kHz / SF8 (any CR).
      Per link (transmitter -> observer) it stores an SNR histogram
      (0.25 dB bins). The transmitter is the last path hop; with an empty path
      it is the originator for ADVERT/ANON_REQ (full pubkey known), otherwise
      the observation is skipped (originator unknown). Per (packet, observer)
      the best SNR over all copies is histogrammed too ("delivery" view).
      Consecutive path hops are counted as repeater adjacencies (no SNR).

  report <links.json.gz> --nodes nodes.json --observers observers.json
         [--center LAT LON] [--radius KM ...] [--shift DB ...]
      Resolves hop prefixes to nodes, computes per-link and per-observer
      statistics and the share of packets that would fall below an SF7
      decode floor. The SF8 floor is taken from the data (low percentile of
      all decoded SNRs) and the SF7 floor = SF8 floor + shift (default 2.5 dB,
      the Semtech SNR-limit step between SF8 and SF7).

Assumption (stated in every report): the physical SNR of a link is the same at
SF7 and SF8 (same 62.5 kHz noise bandwidth); only the demodulation floor moves.
Packets SF8 never decoded are invisible, and a real switch also changes interference,
retries and forwarding, so the shares are threshold EXPOSURE of today's receipts, not a forecast.
"""

import argparse
import gzip
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict

FREQ_MHZ = 869.618
FREQ_TOL_MHZ = 0.0006
PT_ADVERT = 4
PT_ANON_REQ = 7
PT_TRACE = 9
FLOOD_ROUTES = (0, 1)
HOP_RE = re.compile(r"[0-9A-F]{2}|[0-9A-F]{4}|[0-9A-F]{6}")  # TRANSPORT_FLOOD, FLOOD: path = hops so far, last entry = transmitter


def radio_profile(radio, sf_wanted=8):
    """'ok' for 869.618 MHz / 62.5 kHz / wanted SF (any CR); 'other' for a parsed other profile; 'unknown' else."""
    parts = (radio or "").split(",")
    if len(parts) != 4:
        return "unknown"
    try:
        f, bw, sf = float(parts[0]), float(parts[1]), int(parts[2])
    except ValueError:
        return "unknown"
    if abs(f - FREQ_MHZ) <= FREQ_TOL_MHZ and abs(bw - 62.5) < 1e-6 and sf == sf_wanted:
        return "ok"
    return "other"
BIN = 0.25


def snr_bin(snr):
    return round(round(snr / BIN) * BIN, 2)


def originator(record):
    """Full originator pubkey for ADVERT / ANON_REQ, else None."""
    if record.get("payload_type") not in (PT_ADVERT, PT_ANON_REQ):
        return None
    try:
        dec = json.loads(record.get("decoded_json") or "{}")
    except ValueError:
        return None
    key = dec.get("pubKey") if record.get("payload_type") == PT_ADVERT else dec.get("srcPubKey")
    if isinstance(key, str) and re.fullmatch(r"[0-9a-fA-F]{64}", key):
        return key.upper()
    return None


class Acc:
    def __init__(self):
        self.links = defaultdict(Counter)       # (tx, obs) -> Counter(bin)
        self.best = Counter()                   # bin of best SNR per (pkt, obs)
        self.best_obs = defaultdict(Counter)    # obs -> Counter(bin of best SNR per packet)
        self.adj = Counter()                    # (a, b) sorted hop pair -> n
        self.meta = Counter()

    def add_record(self, rec, seen_ids=None):
        """Receipts (best copy per packet+observer) use EVERY eligible SF8 observation; link attribution
        only where the transmitter is known: flood route (not TRACE) with a path -> last hop; empty path
        on ADVERT/ANON_REQ -> originator key. Direct routes with a path carry the REMAINING route -> no link."""
        orig, orig_done = None, False
        per_obs_best = {}
        route = rec.get("route_type")
        ptype = rec.get("payload_type")
        for o in rec.get("observations") or []:
            self.meta["observations"] += 1
            oid = o.get("id")
            if seen_ids is not None and oid is not None:
                if oid in seen_ids:
                    self.meta["dup_observation"] += 1
                    continue
                seen_ids.add(oid)
            if (o.get("direction") or "rx") != "rx":
                self.meta["not_rx"] += 1
                continue
            prof = radio_profile(o.get("radio"))
            if prof != "ok":
                self.meta["radio_" + prof] += 1
                continue
            snr = o.get("snr")
            try:
                snr = float(snr)
            except (TypeError, ValueError):
                self.meta["no_snr"] += 1
                continue
            if not math.isfinite(snr):
                self.meta["no_snr"] += 1
                continue
            obs = (o.get("observer_id") or "").upper()
            if not obs:
                self.meta["no_observer"] += 1
                continue
            b = snr_bin(snr)
            self.meta["receipts_obs"] += 1
            if obs not in per_obs_best or b > per_obs_best[obs]:
                per_obs_best[obs] = b
            try:
                path = json.loads(o.get("path_json") or "[]")
            except ValueError:
                self.meta["bad_path"] += 1
                continue
            if not isinstance(path, list):
                self.meta["bad_path"] += 1
                continue
            path = [str(h).upper() for h in path if h]
            if path and not all(HOP_RE.fullmatch(h) for h in path):
                self.meta["bad_hop"] += 1
                continue
            if path and len({len(h) for h in path}) != 1:
                self.meta["mixed_hop_width"] += 1
                continue
            if path and (route not in FLOOD_ROUTES or ptype == PT_TRACE):
                self.meta["skip_direct_or_trace"] += 1
                continue
            if path:
                for a, c in zip(path, path[1:]):
                    if a != c:
                        self.adj[tuple(sorted((a, c)))] += 1
                tx = path[-1]
            else:
                # empty path: the transmitter is the originator only for ADVERT (never direct-forwarded)
                # or a FLOOD ANON_REQ; an empty DIRECT route is the last forwarded copy (sender unknown)
                if ptype == PT_ANON_REQ and route not in FLOOD_ROUTES:
                    self.meta["skip_unknown_tx"] += 1
                    continue
                if not orig_done:
                    orig = originator(rec)
                    orig_done = True
                if orig is None:
                    self.meta["skip_unknown_tx"] += 1
                    continue
                tx = orig
            self.links[(tx, obs)][b] += 1
            self.meta["used"] += 1
        for o_id, b in per_obs_best.items():
            self.best[b] += 1
            self.best_obs[o_id][b] += 1
            self.meta["deliveries"] += 1

    def merge(self, other):
        for k, c in other.links.items():
            self.links[k].update(c)
        self.best.update(other.best)
        for k, c in other.best_obs.items():
            self.best_obs[k].update(c)
        self.adj.update(other.adj)
        self.meta.update(other.meta)

    def to_json(self):
        return {
            "links": [[tx, obs, {str(b): n for b, n in c.items()}]
                      for (tx, obs), c in self.links.items()],
            "best": {str(b): n for b, n in self.best.items()},
            "best_obs": {o: {str(b): n for b, n in c.items()} for o, c in self.best_obs.items()},
            "adj": [[a, b, n] for (a, b), n in self.adj.items()],
            "meta": dict(self.meta),
        }

    @classmethod
    def from_json(cls, d):
        acc = cls()
        for tx, obs, h in d["links"]:
            acc.links[(tx, obs)] = Counter({float(b): n for b, n in h.items()})
        acc.best = Counter({float(b): n for b, n in d["best"].items()})
        for o, h in d.get("best_obs", {}).items():
            acc.best_obs[o] = Counter({float(b): n for b, n in h.items()})
        acc.adj = Counter({(a, b): n for a, b, n in d["adj"]})
        acc.meta = Counter(d["meta"])
        return acc


def extract_file(path):
    acc = Acc()
    seen = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except ValueError:
                acc.meta["bad_line"] += 1
                continue
            acc.add_record(rec, seen)
    acc.meta["files"] += 1
    return acc


# ---------------------------------------------------------------- statistics

def hist_quantile(hist, q):
    total = sum(hist.values())
    if total == 0:
        return None
    target = q * total
    run = 0
    for b in sorted(hist):
        run += hist[b]
        if run >= target:
            return b
    return max(hist)


def frac_below(hist, thr):
    total = sum(hist.values())
    if total == 0:
        return None
    return sum(n for b, n in hist.items() if b < thr) / total


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def has_pos(n):
    lat, lon = n.get("lat"), n.get("lon")
    return lat is not None and lon is not None and not (abs(lat) < 1e-6 and abs(lon) < 1e-6)


def build_index(nodes, observers):
    """pubkey(upper) -> {name, lat, lon, role}."""
    idx = {}
    for n in nodes:
        pk = (n.get("public_key") or n.get("pubkey") or "").upper()
        if pk:
            idx[pk] = {"name": n.get("name"), "lat": n.get("lat"), "lon": n.get("lon"),
                       "role": n.get("role")}
    for o in observers:
        pk = (o.get("id") or "").upper()
        if not pk:
            continue
        e = idx.setdefault(pk, {"name": o.get("name"), "lat": None, "lon": None, "role": None})
        if not has_pos(e) and o.get("lat") is not None:
            e["lat"], e["lon"] = o.get("lat"), o.get("lon")
        e["observer"] = True
        if not e.get("name"):
            e["name"] = o.get("name")
    return idx


def resolve_prefix(prefix, obs_pk, idx, sorted_keys, cache, strict=False):
    """Resolve a hop prefix seen by an observer to one pubkey or None.

    Unique match wins. Several matches: the one closest to the observer wins
    if it is under half the distance of the runner-up (geo plausibility).
    Repeaters/rooms are preferred over companions when ambiguous.
    """
    if len(prefix) == 64:
        return prefix if prefix in idx else prefix
    key = (prefix, obs_pk, strict)
    if key in cache:
        return cache[key]
    import bisect
    i = bisect.bisect_left(sorted_keys, prefix)
    cands = []
    while i < len(sorted_keys) and sorted_keys[i].startswith(prefix):
        cands.append(sorted_keys[i])
        i += 1
    res = None
    if len(cands) == 1:
        res = cands[0]
    elif len(cands) > 1 and not strict:
        rep = [c for c in cands if (idx[c].get("role") or "").lower() in ("repeater", "room", "room_server")]
        if len(rep) == 1:
            res = rep[0]
        else:
            pool = rep or cands
            o = idx.get(obs_pk)
            if o and has_pos(o) and all(has_pos(idx[c]) for c in pool):
                d = sorted((haversine_km(o["lat"], o["lon"], idx[c]["lat"], idx[c]["lon"]), c)
                           for c in pool if has_pos(idx[c]))
                if len(d) == 1 or (len(d) >= 2 and d[0][0] < 0.5 * d[1][0]):
                    res = d[0][1]
    cache[key] = res
    return res


def analyse(acc, idx, center=None, radii=(), shifts=(2.5,), floor_q=0.005, min_n=5, strict=False):
    sorted_keys = sorted(idx)
    cache = {}
    all_hist = Counter()
    for c in acc.links.values():
        all_hist.update(c)
    if not all_hist:
        raise SystemExit("no SF8 observations with SNR in the input")
    if not 0 <= floor_q <= 1:
        raise SystemExit("floor quantile must be in [0, 1]")
    f8 = hist_quantile(all_hist, floor_q)
    out = {"meta": dict(acc.meta), "sf8_floor": f8, "floor_q": floor_q,
           "snr_all_p": {q: hist_quantile(all_hist, q) for q in (0.01, 0.05, 0.1, 0.25, 0.5)},
           "shifts": {}}

    merged = defaultdict(Counter)
    keyinfo = {}
    res_counts = Counter()
    for (tx, obs), h in acc.links.items():
        n = sum(h.values())
        txpk = resolve_prefix(tx, obs, idx, sorted_keys, cache, strict)
        if txpk is None:
            res_counts["unresolved"] += n
        elif txpk not in idx:
            res_counts["full_key_not_in_catalogue"] += n
        else:
            res_counts["resolved"] += n
        k = (txpk or tx, obs)
        merged[k].update(h)
        keyinfo[k] = (tx, txpk)
    links = []
    for k, h in merged.items():
        tx, txpk = keyinfo[k]
        links.append({"tx": tx, "txpk": txpk, "obs": k[1], "n": sum(h.values()), "hist": h,
                      "p10": hist_quantile(h, 0.1), "p50": hist_quantile(h, 0.5)})
    out["resolution_obs"] = dict(res_counts)
    out["unresolved_obs"] = res_counts["unresolved"]
    out["links_total"] = len(links)

    def dist_center(pk):
        if center is None or pk is None:
            return None
        e = idx.get(pk)
        if not e or not has_pos(e):
            return None
        return haversine_km(center[0], center[1], e["lat"], e["lon"])

    for s in shifts:
        f7 = f8 + s
        res = {"sf7_floor": f7,
               "delivery_loss": frac_below(acc.best, f7),
               "observation_loss": frac_below(all_hist, f7)}
        areas = {}
        for r in list(radii) + [None]:
            sel = [lk for lk in links
                   if r is None or ((dist_center(lk["obs"]) is not None and dist_center(lk["obs"]) <= r))]
            sel_n = [lk for lk in sel if lk["n"] >= min_n]
            h = Counter()
            for lk in sel:
                h.update(lk["hist"])
            obs_best = {}
            for lk in sel_n:
                p50 = lk["p50"]
                if lk["obs"] not in obs_best or p50 > obs_best[lk["obs"]][0]:
                    obs_best[lk["obs"]] = (p50, lk["txpk"] or lk["tx"])
            cls = Counter()
            for lk in sel_n:
                lost = frac_below(lk["hist"], f7)
                cls["lost" if lost >= 0.5 else "degraded" if lost >= 0.1 else "ok"] += 1
            obs_status = Counter()
            for p50, _ in obs_best.values():
                obs_status["cut" if p50 < f7 else "marginal" if p50 < f7 + 3 else "ok"] += 1
            bh = Counter()
            for o_id, c in acc.best_obs.items():
                if r is None or (dist_center(o_id) is not None and dist_center(o_id) <= r):
                    bh.update(c)
            areas["all" if r is None else f"{r:g} km"] = {
                "deliveries": sum(bh.values()),
                "delivery_loss": frac_below(bh, f7) if bh else None,
                "observations": sum(h.values()),
                "observation_loss": frac_below(h, f7) if h else None,
                "links_min_n": len(sel_n),
                "link_classes": dict(cls),
                "observers": len(obs_best),
                "observer_best_link": dict(obs_status),
                "observers_detail": sorted(
                    ((idx.get(o, {}).get("name") or o[:8], p50, idx.get(t, {}).get("name") if t else None)
                     for o, (p50, t) in obs_best.items()), key=lambda x: x[1])[:25],
            }
        res["areas"] = areas
        out["shifts"][str(s)] = res
    out["_links"] = links
    return out


def node_view(out, idx, prefix):
    """Links where the transmitter or the observer matches a pubkey prefix."""
    prefix = prefix.upper()
    rows = []
    for lk in out["_links"]:
        tx = lk["txpk"] or lk["tx"]
        if tx.startswith(prefix) or lk["obs"].startswith(prefix):
            rows.append((idx.get(tx, {}).get("name") or tx[:8],
                         idx.get(lk["obs"], {}).get("name") or lk["obs"][:8],
                         lk["n"], lk["p10"], lk["p50"]))
    return sorted(rows, key=lambda r: -r[2])


def neighbour_samples(lines, since=None, tol_s=150):
    """Unique zero-hop advert SNR samples from the repeater logger 'neighbours' records.

    The CLI table shows the SNR of the LAST advert heard per neighbour; the same advert reappears in
    every poll (10 min) with a growing secs_ago. Samples are de-duplicated globally (all files) by
    (id, snr) with heard times within tol_s; '--since' applies to the heard time. Adverts overwritten
    between two polls are not visible (sampling, not a complete record).
    """
    from datetime import datetime, timedelta
    raw = []
    for line in lines:
        if '"neighbours"' not in line:
            continue
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if r.get("type") != "neighbours":
            continue
        ts = datetime.fromisoformat(r["ts"])
        for n in r.get("neighbours") or []:
            if n.get("snr_db") is None or n.get("secs_ago") is None:
                continue
            heard = ts - timedelta(seconds=n["secs_ago"])
            if since and heard < since:
                continue
            raw.append((n["id"].upper(), heard.timestamp(), float(n["snr_db"])))
    raw.sort()
    per = defaultdict(list)
    last = {}
    for nid, t, snr in raw:
        prev = last.get((nid, snr))
        if prev is not None and t - prev <= tol_s:
            last[(nid, snr)] = t
            continue
        last[(nid, snr)] = t
        per[nid].append(snr)
    return per


def cmd_neigh(args):
    from datetime import datetime
    since = datetime.fromisoformat(args.since) if args.since else None
    lines = []
    for p in args.logs:
        with open(p) as fh:
            lines.extend(fh.readlines())
    per = neighbour_samples(lines, since)
    idx = {}
    if args.nodes:
        nodes = json.load(open(args.nodes))
        nodes = nodes.get("nodes", nodes) if isinstance(nodes, dict) else nodes
        idx = {(n.get("public_key") or "").upper(): n for n in nodes}
    keys = sorted(idx)
    print("1863 neighbour table: unique zero-hop advert SNR samples (as heard BY 1863)")
    allv = [x for v in per.values() for x in v]
    for thr in args.thr:
        nb = sum(1 for v in per.values() if statistics.median(v) < thr)
        print(f"threshold {thr} dB: samples below {sum(1 for x in allv if x < thr)}/{len(allv)}, "
              f"neighbours with median below {nb}/{len(per)}")
    rows = []
    for nid, v in per.items():
        v = sorted(v)
        name, dist = "", None
        cands = [k for k in keys if k.startswith(nid)]
        if len(cands) == 1:
            n = idx[cands[0]]
            name = n.get("name") or ""
            if args.center and has_pos(n):
                dist = haversine_km(args.center[0], args.center[1], n["lat"], n["lon"])
        rows.append((statistics.median(v), nid, name, len(v), v[0], dist))
    for med, nid, name, cnt, mn, dist in sorted(rows):
        print(f"  {nid} {name[:28]:28s} n={cnt:3d} min={mn:6.2f} median={med:6.2f} dist={fmt(dist)} km")


# ------------------------------------------------------------------- CLI

def cmd_extract(args):
    acc = Acc()
    if args.workers > 1 and len(args.dumps) > 1:
        from multiprocessing import Pool
        with Pool(args.workers) as pool:
            for part in pool.imap_unordered(extract_file, args.dumps):
                acc.merge(part)
                print(f"merged {sum(part.links[k].total() for k in part.links)} obs", file=sys.stderr)
    else:
        for p in args.dumps:
            acc.merge(extract_file(p))
            print(f"done {p}", file=sys.stderr)
    with gzip.open(args.out, "wt") as fh:
        json.dump(acc.to_json(), fh)
    print(json.dumps(dict(acc.meta)))


def fmt(v, nd=1):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def cmd_report(args):
    with gzip.open(args.links, "rt") as fh:
        acc = Acc.from_json(json.load(fh))
    nodes = json.load(open(args.nodes))
    nodes = nodes.get("nodes", nodes) if isinstance(nodes, dict) else nodes
    observers = json.load(open(args.observers))
    observers = observers.get("observers", observers) if isinstance(observers, dict) else observers
    idx = build_index(nodes, observers)
    center = tuple(args.center) if args.center else None
    if args.radius and center is None:
        raise SystemExit("--radius needs --center")
    out = analyse(acc, idx, center, args.radius, args.shift, args.floor_q, args.min_n)
    print("Assumption: same link SNR at SF7 and SF8; SF7 scenario threshold = 0.5th percentile of decoded SF8 copies + shift. "
          "Shares = THRESHOLD EXPOSURE of today's receipts (frozen hard-threshold model), not a delivery forecast.")
    print(f"meta: {out['meta']}")
    print(f"SF8 decoded-copy percentile (q={out['floor_q']}): {out['sf8_floor']} dB; SNR quantiles attributed copies: {out['snr_all_p']}")
    print(f"links: {out['links_total']}, attributed copies with unresolved transmitter (heuristic rule): {out['unresolved_obs']}")
    for s, res in out["shifts"].items():
        print(f"\n=== shift {s} dB -> SF7 scenario threshold {res['sf7_floor']} dB ===")
        print(f"deliveries (packet,observer; best copy) below threshold: {fmt(100*res['delivery_loss'])} %")
        print(f"observations below threshold: {fmt(100*res['observation_loss'])} %")
        for area, a in res["areas"].items():
            print(f"-- {area}: deliveries {a['deliveries']} below threshold {fmt(100*a['delivery_loss']) if a['delivery_loss'] is not None else '-'} %, obs {a['observations']}, below threshold {fmt(100*a['observation_loss']) if a['observation_loss'] is not None else '-'} %, "
                  f"links(n>={args.min_n}) {a['links_min_n']} {a['link_classes']}, observers {a['observers']} best-link {a['observer_best_link']}")
            if args.details:
                for name, p50, via in a["observers_detail"]:
                    print(f"     {name}: best link median {p50} dB via {via}")
    for pfx in args.node or []:
        print(f"\n=== links involving {pfx} (tx, observer, n, p10, p50) ===")
        for r in node_view(out, idx, pfx)[:40]:
            print("  ", r)
    if args.json_out:
        dump = {k: v for k, v in out.items() if k != "_links"}
        dump["links"] = [{"tx": lk["txpk"] or lk["tx"], "obs": lk["obs"], "n": lk["n"],
                          "p10": lk["p10"], "p50": lk["p50"],
                          "hist": {str(b): c for b, c in lk["hist"].items()}} for lk in out["_links"]]
        with open(args.json_out, "w") as fh:
            json.dump(dump, fh)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("dumps", nargs="+")
    e.add_argument("--out", required=True)
    e.add_argument("--workers", type=int, default=1)
    r = sub.add_parser("report")
    r.add_argument("links")
    r.add_argument("--nodes", required=True)
    r.add_argument("--observers", required=True)
    r.add_argument("--center", type=float, nargs=2)
    r.add_argument("--radius", type=float, nargs="*", default=[])
    r.add_argument("--shift", type=float, nargs="*", default=[2.5])
    r.add_argument("--floor-q", type=float, default=0.005)
    r.add_argument("--min-n", type=int, default=5)
    r.add_argument("--node", nargs="*")
    r.add_argument("--details", action="store_true")
    r.add_argument("--json-out")
    g = sub.add_parser("neigh")
    g.add_argument("logs", nargs="+")
    g.add_argument("--since")
    g.add_argument("--nodes")
    g.add_argument("--center", type=float, nargs=2)
    g.add_argument("--thr", type=float, nargs="*", default=[-7.5, -10.0])
    args = ap.parse_args(argv)
    if args.cmd == "extract":
        cmd_extract(args)
    elif args.cmd == "neigh":
        cmd_neigh(args)
    else:
        cmd_report(args)


if __name__ == "__main__":
    main()
