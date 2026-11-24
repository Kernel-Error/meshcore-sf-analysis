#!/usr/bin/env python3
"""SF7/SF8 assessment for the Rheinland region: metrics + figures (figure text German).

Usage:
  sf7_graphics.py <outdir> --links links.json.gz --nodes nodes.json --observers observers.json
                  --raklogs rak-*.jsonl [--neigh-since ISO] [--load-since ISO --load-until ISO]

Writes sf7-metrics.json and:
  sf7-01-tradeoff.png        modelled link budget vs airtime (100 B PHY payload) per SF/CR/BW
  sf7-02-snr-cdf.png         share of decoded SF8 copies strictly below x (observers <= 60 km; 1863 neighbour adverts)
  sf7-03-exposure-vs-threshold.png  share of today's receipts whose best copy is below an assumed SF7 threshold
  sf7-04-observer-coverage.png      where the observers are (who measures what), 15/60 km circles, 1863
  sf7-05-distance-classes.png       links by distance: share of copies below the SF7 scenario threshold
  sf7-06-channel-model.png   pure-ALOHA illustration with ASSUMED offered load (no measured anchor)
  sf7-07-neigh-distance.png  1863 neighbour advert SNR (median) vs distance, singletons marked

Terms: 'SF8 percentile' = 0.5th percentile of decoded SF8 copies (descriptive, not a demodulation limit).
'SF7 scenario' = that percentile + 2.5 dB (nominal Semtech step; +3 dB sensitivity); -7.5 dB = datasheet reference.
Shares are THRESHOLD EXPOSURE of today's receipts in a frozen hard-threshold model, not a delivery forecast:
packets SF8 never decoded are invisible, and a real switch also changes interference, retries and forwarding.
Equal link SNR at SF7/SF8 and equal CR sensitivity are modelling assumptions. ALOHA ignores LBT, capture,
hidden nodes, retries and correlated flood starts; it is a relative illustration only.
"""
import argparse
import gzip
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lora_airtime as la  # noqa: E402
import sf7_link_analysis as sla  # noqa: E402

POS_1863 = (50.577972, 6.949222)  # public offset position of repeater 1863
CENTER = POS_1863                  # analysis centre (separate constant on purpose)
OKABE = {"blue": "#0072B2", "orange": "#E69F00", "vermilion": "#D55E00", "green": "#009E73",
         "sky": "#56B4E9", "grey": "#777777", "purple": "#CC79A7"}
MIN_PT = 21
SF_STEP_DB = 2.5
PAYLOAD = 100
DIST_BINS = [(0, 2), (2, 10), (10, 20), (20, 40), (40, 80)]
LINK_MIN_N = 20

# (label, sf, bw_hz, cr 1..4, note) note: None = nominal fit at 869.618 / 62.5 kHz
CONFIGS = [
    ("SF6/CR5", 6, 62500, 1, None),
    ("SF7/CR5", 7, 62500, 1, None),
    ("SF7/CR8", 7, 62500, 4, None),
    ("SF8/CR5", 8, 62500, 1, None),
    ("SF8/CR8", 8, 62500, 4, None),
    ("SF9/CR5", 9, 62500, 1, None),
    ("SF9/CR8", 9, 62500, 4, None),
    ("SF10/CR5", 10, 62500, 1, None),
    ("SF8/CR5 125k", 8, 125000, 1, "band"),     # 869.618 +/- 62.5 kHz exceeds 869.65 MHz
    ("SF7/CR5 31k", 7, 31250, 1, "clock"),      # +/-7.8 kHz offset rule ~ +/-4.5 ppm per node
]
REF = ("SF8/CR8", 8, 62500, 4)


def de_num(v, nd=1):
    """German number format with a real minus sign: -7.5 -> '−7,5'."""
    return f"{v:.{nd}f}".replace("-", "−").replace(".", ",")


def config_table(payload=PAYLOAD):
    ref_air = la.airtime_s(payload, REF[1], REF[2], REF[3])
    ref_sens = la.sensitivity_dbm(REF[1], REF[2])
    rows = []
    for label, sf, bw, cr, note in CONFIGS:
        air = la.airtime_s(payload, sf, bw, cr)
        sens = la.sensitivity_dbm(sf, bw)
        rows.append({"label": label, "sf": sf, "bw": bw, "cr": cr + 4, "note": note,
                     "airtime_ms": round(air * 1000, 1), "airtime_rel": round(air / ref_air, 3),
                     "snr_limit_db": la.SNR_LIMIT_DB[sf], "sensitivity_dbm": round(sens, 1),
                     "budget_vs_ref_db": round(ref_sens - sens, 1)})
    return rows


def aloha_success(load):
    """Pure ALOHA: probability of zero overlap for independent starts and equal frame length, offered load G."""
    return math.exp(-2.0 * load)


def channel_load(lines, since=None, until=None, min_hours=23.5):
    """Per UTC day: decoded RX airtime share and own TX share from logger stats.

    Records are sorted by timestamp and de-duplicated; deltas only between consecutive samples of the same boot
    (uptime increases), at most 1 h apart and consistent with wall time; days with < min_hours coverage dropped.
    A delta is assigned to the day of its end sample (midnight split error < 1 poll interval).
    """
    from datetime import datetime
    recs = {}
    for line in lines:
        if '"stats"' not in line:
            continue
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if r.get("type") != "stats":
            continue
        core, radio = r.get("core") or {}, r.get("radio") or {}
        up, rx, tx = core.get("uptime_secs"), radio.get("rx_air_secs"), radio.get("tx_air_secs")
        if up is None or rx is None or tx is None:
            continue
        ts = datetime.fromisoformat(r["ts"])
        if (since and ts < since) or (until and ts >= until):
            continue
        recs[ts] = (up, rx, tx)
    acc = defaultdict(lambda: [0.0, 0.0, 0.0])
    prev = None
    for ts in sorted(recs):
        up, rx, tx = recs[ts]
        if prev:
            pts, (pup, prx, ptx) = prev
            dup, dwall = up - pup, (ts - pts).total_seconds()
            if 0 < dup < 3600 and rx >= prx and tx >= ptx and abs(dup - dwall) < 120:
                a = acc[ts.date().isoformat()]
                a[0] += dup
                a[1] += rx - prx
                a[2] += tx - ptx
        prev = (ts, (up, rx, tx))
    return {d: {"rx_share": v[1] / v[0], "tx_share": v[2] / v[0], "hours": v[0] / 3600}
            for d, v in sorted(acc.items()) if v[0] >= min_hours * 3600}


def loss_curve(hist, thresholds):
    return [(f, sla.frac_below(hist, f)) for f in thresholds]


def in_radius(idx, pk, center, radius):
    if radius is None:
        return True
    e = idx.get(pk)
    return bool(e and sla.has_pos(e) and sla.haversine_km(center[0], center[1], e["lat"], e["lon"]) <= radius)


def area_hists(acc, idx, center, radius):
    """(observation hist, best-copy receipt hist, observers used) for observers within radius km (None = all)."""
    obs_h, del_h, used = Counter(), Counter(), set()
    for (tx, o), h in acc.links.items():
        if in_radius(idx, o, center, radius):
            obs_h.update(h)
    for o, h in acc.best_obs.items():
        if in_radius(idx, o, center, radius):
            del_h.update(h)
            used.add(o)
    return obs_h, del_h, used


def equal_observer_exposure(acc, observers, thr):
    vals = {o: sla.frac_below(acc.best_obs[o], thr) for o in observers if acc.best_obs.get(o)}
    vals = {o: v for o, v in vals.items() if v is not None}
    if not vals:
        return None
    return {"mean": statistics.mean(vals.values()), "min": min(vals.values()), "max": max(vals.values()),
            "n": len(vals)}


def link_class(share):
    return "lost" if share >= 0.5 else "degraded" if share >= 0.1 else "ok"


def distance_classes(link_rows, bins=DIST_BINS, min_n=LINK_MIN_N):
    """Per distance bin: link counts per class (<10 %, 10 to <50 %, >=50 % of copies below threshold)."""
    out = []
    for lo, hi in bins:
        sel = [r for r in link_rows if lo <= r["km"] < hi and r["n"] >= min_n]
        c = Counter(link_class(r["loss"]) for r in sel)
        n_obs = sum(r["n"] for r in sel)
        out.append({"bin": f"{lo}–{hi} km", "links": len(sel), "ok": c["ok"], "degraded": c["degraded"],
                    "lost": c["lost"], "obs_share": sum(r["n"] * r["loss"] for r in sel) / n_obs if n_obs else None})
    return out


def share_ge10(link_rows, lo=2, hi=80, min_n=LINK_MIN_N):
    sel = [r for r in link_rows if lo <= r["km"] < hi and r["n"] >= min_n]
    k = sum(1 for r in sel if r["loss"] >= 0.1)
    return {"links": len(sel), "ge10": k, "share": k / len(sel) if sel else None}


def cdf_below(hist, xs):
    tot = sum(hist.values())
    return [100 * sum(n for b, n in hist.items() if b < x) / tot for x in xs] if tot else []


# ------------------------------------------------------------------ figures

def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": MIN_PT, "axes.titlesize": MIN_PT + 2, "axes.labelsize": MIN_PT,
                         "xtick.labelsize": MIN_PT - 2, "ytick.labelsize": MIN_PT - 2, "legend.fontsize": MIN_PT - 3,
                         "figure.dpi": 100})
    return plt


def de_axis(plt, axis, nd=0):
    axis.set_major_formatter(plt.FuncFormatter(lambda v, _: de_num(v, nd)))


def fig_tradeoff(plt, rows, path):
    fig, ax = plt.subplots(figsize=(12, 9))
    for r in rows:
        col = OKABE["vermilion"] if r["label"] == "SF8/CR8" else OKABE["blue"] if r["sf"] == 7 else \
            OKABE["green"] if r["sf"] == 9 else OKABE["grey"]
        hollow = r["note"] is not None
        ax.scatter(r["airtime_ms"], r["budget_vs_ref_db"], s=260, color="white" if hollow else col,
                   edgecolors=col, linewidths=3, zorder=3)
        lab = r["label"] + (" (Bandgrenze)" if r["note"] == "band" else " (Quarz)" if r["note"] == "clock" else "")
        ax.annotate(lab, (r["airtime_ms"], r["budget_vs_ref_db"]), textcoords="offset points",
                    xytext=(10, -26 if r["note"] == "band" or r["label"] == "SF8/CR5" else 8), fontsize=MIN_PT - 4)
    ax.axhline(0, color=OKABE["grey"], lw=1, ls=":")
    ax.set_xscale("log")
    ticks = [200, 400, 800, 1600, 3200]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,}".replace(",", ".") for t in ticks])
    ax.minorticks_off()
    de_axis(plt, ax.yaxis)
    ax.set_xlabel(f"Sendezeit für {PAYLOAD} Byte PHY-Nutzlast [ms]")
    ax.set_ylabel("Modelliertes Linkbudget ggü. SF8/CR8 [dB]")
    ax.set_title("Sendezeit gegen Reichweitenreserve")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_snr_cdf(plt, obs_h, neigh_vals, lines, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    xs = [x / 4 for x in range(-64, 1)]
    ax.step(xs, cdf_below(obs_h, xs), where="post", color=OKABE["blue"], lw=3, label="Beobachter ≤ 60 km")
    if neigh_vals:
        ax.step(xs, cdf_below(Counter(neigh_vals), xs), where="post", color=OKABE["orange"], lw=3,
                label="1863 hört Nachbarn")
    for val, lab, c, ls in lines:
        ax.axvline(val, color=c, ls=ls, lw=2.5, label=lab)
    ax.set_xlim(-16, 0)
    ax.set_ylim(0, 40)
    de_axis(plt, ax.xaxis)
    ax.set_xlabel("SNR [dB]")
    ax.set_ylabel("Kopien unter diesem SNR [%]")
    ax.set_title("Wie viel liegt nahe der Grenze?")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=MIN_PT - 5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_exposure(plt, curves, f8, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    cols = [OKABE["blue"], OKABE["orange"], OKABE["green"]]
    for (lab, pts), c in zip(curves.items(), cols):
        ax.step([p[0] for p in pts], [100 * p[1] for p in pts], where="post", color=c, lw=3, label=lab)
    ax.axvspan(f8 + 2.5, f8 + 3.0, color=OKABE["vermilion"], alpha=0.15, label="SF7-Szenario (+2,5…3 dB)")
    ax.axvline(-7.5, color=OKABE["vermilion"], ls=":", lw=2.5, label="SF7 Datenblatt (−7,5 dB)")
    de_axis(plt, ax.xaxis)
    de_axis(plt, ax.yaxis)
    ax.set_xlabel("angenommene SF7-Schwelle [dB SNR]")
    ax.set_ylabel("Empfänge, beste Kopie darunter [%]")
    ax.set_title("Heutige Empfänge unter einer SF7-Schwelle")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=MIN_PT - 5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_coverage(plt, obs_pts, rep_pts, path, radius_km=60):
    fig, ax = plt.subplots(figsize=(12, 12))
    kx = math.cos(math.radians(CENTER[0]))
    ax.scatter([p[1] for p in rep_pts], [p[0] for p in rep_pts], s=14, color=OKABE["grey"], alpha=0.6,
               label="Repeater (Katalog)")
    ax.scatter([p[1] for p in obs_pts], [p[0] for p in obs_pts], s=220, marker="^", color=OKABE["orange"],
               edgecolors="black", zorder=4, label=f"Beobachter mit Position ({len(obs_pts)})")
    for r, ls in ((15, "--"), (60, "-")):
        t = [i * 2 * math.pi / 180 for i in range(181)]
        ax.plot([CENTER[1] + r / 111.0 / kx * math.cos(a) for a in t], [CENTER[0] + r / 111.0 * math.sin(a) for a in t],
                color=OKABE["blue"], ls=ls, lw=2)
        ax.text(CENTER[1], CENTER[0] + r / 111.0 + 0.01, f"{r} km", color=OKABE["blue"], ha="center",
                fontsize=MIN_PT - 4)
    ax.scatter([POS_1863[1]], [POS_1863[0]], marker="*", s=700, color=OKABE["blue"], zorder=5, label="1863")
    d = (radius_km + 8) / 111.0
    ax.set_ylim(CENTER[0] - d, CENTER[0] + d)
    ax.set_xlim(CENTER[1] - d / kx, CENTER[1] + d / kx)
    ax.set_aspect(1 / kx)
    de_axis(plt, ax.xaxis, 1)
    de_axis(plt, ax.yaxis, 1)
    ax.set_xlabel("Länge [°]")
    ax.set_ylabel("Breite [°]")
    ax.set_title("Wo gemessen wird")
    ax.legend(loc="lower left", fontsize=MIN_PT - 5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_distance(plt, rows, thr, path):
    fig, ax = plt.subplots(figsize=(12, 9))
    x = list(range(len(rows)))
    tot = [max(1, r["links"]) for r in rows]
    ok = [100 * r["ok"] / t for r, t in zip(rows, tot)]
    de = [100 * r["degraded"] / t for r, t in zip(rows, tot)]
    lo = [100 * r["lost"] / t for r, t in zip(rows, tot)]
    ax.bar(x, ok, color=OKABE["green"], label="< 10 %")
    ax.bar(x, de, bottom=ok, color=OKABE["orange"], label="10 bis < 50 %")
    ax.bar(x, lo, bottom=[a + b for a, b in zip(ok, de)], color=OKABE["vermilion"], label="≥ 50 %")
    for i, r in enumerate(rows):
        ax.text(i, 101, f"n={r['links']}", ha="center", va="bottom", fontsize=MIN_PT - 5)
    ax.set_xticks(x)
    ax.set_xticklabels([r["bin"].replace(" km", "") for r in rows], fontsize=MIN_PT - 3)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Anteil Verbindungen [%]")
    ax.set_xlabel("Entfernung Sender → Beobachter [km]")
    ax.set_title(f"Kopien unter {de_num(thr)} dB je Verbindung")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=MIN_PT - 4, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_channel(plt, rows, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    ref = next(r for r in rows if r["label"] == "SF8/CR8")
    xs = [i / 100 for i in range(5, 61)]
    for lab, col, ls in (("SF9/CR8", OKABE["green"], "--"), ("SF8/CR8", OKABE["vermilion"], "-"),
                         ("SF8/CR5", OKABE["orange"], "-"), ("SF7/CR8", OKABE["sky"], "-"),
                         ("SF7/CR5", OKABE["blue"], "-")):
        r = next(x for x in rows if x["label"] == lab)
        k = r["airtime_rel"] / ref["airtime_rel"]
        ax.plot([100 * g for g in xs], [100 * aloha_success(g * k) for g in xs], color=col, ls=ls, lw=3, label=lab)
    ax.set_xlabel("angenommene Last bei SF8/CR8 [%]")
    ax.set_ylabel("ohne Überlappung [%] (Modell)")
    ax.set_title("Gleicher Verkehr, andere Sendezeit")
    ax.set_xlim(5, 60)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=MIN_PT - 4)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_neigh_distance(plt, pts, lines, n_total, path):
    fig, ax = plt.subplots(figsize=(12, 8))
    multi = [p for p in pts if p[2] > 1]
    single = [p for p in pts if p[2] == 1]
    ax.scatter([p[0] for p in multi], [p[1] for p in multi], s=[30 + 3 * min(p[2], 60) for p in multi],
               color=OKABE["blue"], alpha=0.7, label="Median (Größe ~ Anzahl, max. 60)")
    ax.scatter([p[0] for p in single], [p[1] for p in single], s=90, facecolors="white",
               edgecolors=OKABE["blue"], linewidths=2, label="nur 1 Messwert")
    for val, lab, c, ls in lines:
        ax.axhline(val, color=c, ls=ls, lw=2.5, label=lab)
    de_axis(plt, ax.yaxis)
    ax.set_xlabel("Entfernung zu 1863 [km]")
    ax.set_ylabel("SNR bei 1863 [dB]")
    ax.set_title(f"Wen 1863 wie gut hört ({len(pts)} von {n_total} mit Position)")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=MIN_PT - 6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------------ main

def load_list(p, key):
    d = json.load(open(p))
    return d.get(key, d) if isinstance(d, dict) else d


def link_rows_for(out, idx, thr):
    rows = []
    for lk in out["_links"]:
        o, t = idx.get(lk["obs"]), idx.get(lk["txpk"]) if lk["txpk"] else None
        if o and t and sla.has_pos(o) and sla.has_pos(t):
            rows.append({"km": sla.haversine_km(o["lat"], o["lon"], t["lat"], t["lon"]), "n": lk["n"],
                         "loss": sla.frac_below(lk["hist"], thr), "obs": lk["obs"]})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--links", required=True)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--observers", required=True)
    ap.add_argument("--raklogs", nargs="+", required=True)
    ap.add_argument("--neigh-since", default="2026-10-27T17:36:39+00:00")
    ap.add_argument("--load-since", default="2026-11-12T00:00:00+00:00")
    ap.add_argument("--load-until", default="2026-11-24T00:00:00+00:00")
    ap.add_argument("--floor-q", type=float, default=0.005)
    args = ap.parse_args(argv)
    from datetime import datetime
    os.makedirs(args.outdir, exist_ok=True)

    with gzip.open(args.links, "rt") as fh:
        acc = sla.Acc.from_json(json.load(fh))
    idx = sla.build_index(load_list(args.nodes, "nodes"), load_list(args.observers, "observers"))
    all_h = Counter()
    for h in acc.links.values():
        all_h.update(h)
    if not all_h:
        raise SystemExit("no attributed SF8 copies in the input")
    f8 = sla.hist_quantile(all_h, args.floor_q)
    f7, f7b = f8 + SF_STEP_DB, f8 + 3.0

    # 1863: neighbours (one global de-dup pass) + channel load
    lines = []
    for p in args.raklogs:
        with open(p) as fh:
            lines.extend(fh.readlines())
    neigh = sla.neighbour_samples(lines, datetime.fromisoformat(args.neigh_since))
    loads = channel_load(lines, datetime.fromisoformat(args.load_since), datetime.fromisoformat(args.load_until))

    rows = config_table()
    metrics = {"terms": __doc__.split("Terms:")[1].strip(), "sf8_p0_5_db": f8, "floor_quantile": args.floor_q,
               "thresholds": {"sf7_scenario": f7, "sf7_plus3": f7b, "datasheet": -7.5},
               "snr_quantiles_attributed": {q: sla.hist_quantile(all_h, q) for q in (0.005, 0.01, 0.05, 0.1, 0.5)},
               "configs": rows, "meta": dict(acc.meta)}
    if loads:
        rx = [v["rx_share"] for v in loads.values()]
        tx = [v["tx_share"] for v in loads.values()]
        metrics["decoded_rx_1863"] = {"days": len(loads), "min": min(rx), "median": statistics.median(rx),
                                      "max": max(rx), "tx_min": min(tx), "tx_max": max(tx),
                                      "note": "decoded RX airtime estimated by the firmware with its own CR5; "
                                              "not RF occupancy, not offered load", "per_day": loads}

    # receipts (best copy) per area
    radii = {"Todenfeld ≤ 15 km": 15, "Rheinland ≤ 60 km": 60, "alle Beobachter": None}
    grid = [round(f8 + i * 0.25, 2) for i in range(0, 4 * 8 + 1)]
    curves, area_m = {}, {}
    obs60 = Counter()
    for lab, r in radii.items():
        oh, dh, used = area_hists(acc, idx, CENTER, r)
        if r == 60:
            obs60 = oh
        if not dh:
            area_m[lab] = {"observers": 0, "status": "no data"}
            continue
        curves[f"{lab} ({len(used)})"] = loss_curve(dh, grid)
        area_m[lab] = {"observers": len(used),
                       "observer_names": sorted(idx.get(o, {}).get("name") or o[:8] for o in used) if r == 15 else None,
                       "copies": sum(oh.values()), "receipts": sum(dh.values()),
                       "receipt_exposure": {"sf7_scenario": sla.frac_below(dh, f7), "sf7_plus3": sla.frac_below(dh, f7b),
                                            "datasheet": sla.frac_below(dh, -7.5)},
                       "copy_exposure": {"sf7_scenario": sla.frac_below(oh, f7), "datasheet": sla.frac_below(oh, -7.5)},
                       "equal_observer_sf7_scenario": equal_observer_exposure(acc, used, f7)}
    metrics["areas"] = area_m
    metrics["observers_without_position"] = sum(1 for o in acc.best_obs if not (idx.get(o) and sla.has_pos(idx[o])))
    metrics["observers_total"] = len(acc.best_obs)

    # links: strict (catalogue-unique) primary, heuristic as sensitivity
    out_s = sla.analyse(acc, idx, CENTER, [], [SF_STEP_DB], args.floor_q, strict=True)
    out_h = sla.analyse(acc, idx, CENTER, [], [SF_STEP_DB], args.floor_q, strict=False)
    lr_s = link_rows_for(out_s, idx, f7)
    lr_h = link_rows_for(out_h, idx, f7)
    lr60 = [r for r in lr_s if in_radius(idx, r["obs"], CENTER, 60)]
    metrics["resolution_strict"] = out_s["resolution_obs"]
    metrics["resolution_heuristic"] = out_h["resolution_obs"]
    metrics["links"] = {
        "rule": f"directed transmitter->observer, both positions known, n >= {LINK_MIN_N} attributed copies",
        "strict_all": share_ge10(lr_s), "strict_rx60": share_ge10(lr60),
        "heuristic_all": share_ge10(lr_h),
        "strict_all_plus3": share_ge10(link_rows_for(out_s, idx, f7b)),
        "strict_all_datasheet": share_ge10(link_rows_for(out_s, idx, -7.5)),
        "distance_classes_strict": distance_classes(lr_s),
        "distance_classes_strict_datasheet": distance_classes(link_rows_for(out_s, idx, -7.5)),
        "beyond_80km": sum(1 for r in lr_s if r["km"] >= 80 and r["n"] >= LINK_MIN_N)}

    # 1863 neighbours
    keys = sorted(idx)
    pts, nvals = [], []
    for nid, v in neigh.items():
        v = sorted(v)
        nvals.extend(v)
        c = [k for k in keys if k.startswith(nid)]
        if len(c) == 1 and sla.has_pos(idx[c[0]]):
            e = idx[c[0]]
            pts.append((sla.haversine_km(POS_1863[0], POS_1863[1], e["lat"], e["lon"]), statistics.median(v), len(v)))
    nh = Counter(nvals)
    metrics["neigh_1863"] = {
        "neighbours": len(neigh), "samples": len(nvals), "with_position": len(pts),
        "singletons": sum(1 for v in neigh.values() if len(v) == 1),
        "samples_below": {"sf7_scenario": sla.frac_below(nh, f7), "datasheet": sla.frac_below(nh, -7.5)},
        "medians_below": {"sf7_scenario": sum(1 for v in neigh.values() if statistics.median(v) < f7),
                          "datasheet": sum(1 for v in neigh.values() if statistics.median(v) < -7.5)},
        "dedup": "global over all logs; same (id, snr) within 150 s = same advert; cutoff on heard time"}

    obs_pts = [(idx[o]["lat"], idx[o]["lon"]) for o in acc.best_obs if idx.get(o) and sla.has_pos(idx[o])]
    rep_pts = sorted({(e["lat"], e["lon"]) for e in idx.values() if sla.has_pos(e) and (e.get("role") or "") == "repeater"})
    metrics["coverage"] = {"observers_plotted": len(obs_pts), "repeater_coordinates": len(rep_pts)}

    with open(os.path.join(args.outdir, "sf7-metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1, ensure_ascii=False, default=str)

    lines_cdf = [(f8, f"0,5-Perzentil SF8 ({de_num(f8)} dB)", OKABE["grey"], "-"),
                 (f7, f"SF7-Szenario ({de_num(f7)} dB)", OKABE["vermilion"], "--"),
                 (-7.5, "SF7 Datenblatt (−7,5 dB)", OKABE["vermilion"], ":")]
    plt = setup_mpl()
    fig_tradeoff(plt, rows, os.path.join(args.outdir, "sf7-01-tradeoff.png"))
    fig_snr_cdf(plt, obs60, nvals, lines_cdf, os.path.join(args.outdir, "sf7-02-snr-cdf.png"))
    if curves:
        fig_exposure(plt, curves, f8, os.path.join(args.outdir, "sf7-03-exposure-vs-threshold.png"))
    fig_coverage(plt, obs_pts, rep_pts, os.path.join(args.outdir, "sf7-04-observer-coverage.png"))
    fig_distance(plt, metrics["links"]["distance_classes_strict"], f7,
                 os.path.join(args.outdir, "sf7-05-distance-classes.png"))
    fig_channel(plt, rows, os.path.join(args.outdir, "sf7-06-channel-model.png"))
    fig_neigh_distance(plt, pts, lines_cdf[1:], len(neigh), os.path.join(args.outdir, "sf7-07-neigh-distance.png"))
    show = {k: metrics[k] for k in ("sf8_p0_5_db", "thresholds", "areas", "links", "neigh_1863", "coverage",
                                    "resolution_strict", "observers_without_position")}
    show["decoded_rx_1863"] = {k: v for k, v in metrics.get("decoded_rx_1863", {}).items() if k != "per_day"}
    print(json.dumps(show, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
