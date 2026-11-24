# meshcore-sf-analysis

Small Python tools to estimate what a MeshCore region switch from **SF8 to SF7** (869.618 MHz, 62.5 kHz) would mean, based on the link SNRs a [CoreScope](https://github.com/Kpa-clawbot/CoreScope) instance has observed. It was written for the Rheinland region in Germany, but it works with any CoreScope `packets-expanded` day archive.

Background and the results for Rheinland are in the blog post (German): [FILL: link to the blog post]

## What it does

- `sf7_link_analysis.py extract`: streams CoreScope day dumps (stdlib only, runs fine on a Raspberry Pi).
  - Keeps the RX observations on 869.618 MHz / 62.5 kHz / SF8 (any CR).
  - Builds an SNR histogram for every link (transmitter -> observer), plus the best copy for every (packet, observer).
  - The transmitter of a copy is the last hop of a flood path. An empty path counts as the originator only for ADVERT and flood ANON_REQ. Direct routes and TRACE are not attributed.
- `sf7_link_analysis.py report`: resolves hop prefixes to nodes and prints the threshold exposure per area.
- `sf7_link_analysis.py neigh`: reads the `neighbours` records of a repeater logger and computes per-neighbour advert SNR statistics, de-duplicated across polls.
- `sf7_graphics.py`: writes the metrics JSON and the figures (German labels):
  - airtime vs link budget
  - SNR distribution
  - exposure vs assumed threshold
  - observer coverage
  - links by distance
  - ALOHA illustration
  - neighbour SNR vs distance
- `lora_airtime.py`: time on air (Semtech sx126x driver formula, MeshCore preamble: 32 symbols at SF<=8, 16 above) and a nominal sensitivity estimate.

`results/sf7-metrics-rheinland-2026-11.json` is the metrics file of the Rheinland run: 18 days in November 2026, about 14.9 M observations from 65 observers.

## Usage

```sh
# 1) extract (on the machine that holds the CoreScope archives)
python3 sf7_link_analysis.py extract packets-expanded-2026-11-*.jsonl.gz --out links.json.gz --workers 3

# 2) node + observer snapshot from the CoreScope API (some endpoints want a User-Agent)
curl -A curl/8.5 'https://<corescope>/api/nodes?limit=5000' -o nodes.json
curl -A curl/8.5 'https://<corescope>/api/observers' -o observers.json

# 3) report around a centre point
python3 sf7_link_analysis.py report links.json.gz --nodes nodes.json --observers observers.json \
    --center 50.578 6.949 --radius 15 60 --shift 2.5 3

# 4) figures + metrics (needs matplotlib); --raklogs = repeater logger JSONL with stats/neighbours records
python3 sf7_graphics.py out/ --links links.json.gz --nodes nodes.json --observers observers.json --raklogs rak-*.jsonl
```

`sf7_graphics.py` has its centre point, date windows and the repeater used for the neighbour figure hard-coded for the Rheinland run (`POS_1863`, the `--*-since` defaults). Change them for your own region.

Tests: `python3 -m unittest` (no matplotlib needed for the tests)

## Assumptions and limits (please read before quoting numbers)

- **The link SNR is assumed to be the same at SF7 and SF8**, since the noise bandwidth is the same. Only the demodulation threshold moves. Equal sensitivity for CR5 and CR8 is also assumed; the coding gain is not modelled.
- **Threshold**: the 0.5th percentile of the decoded SF8 copies is descriptive, not a demodulation limit. The SF7 scenario = that percentile + 2.5 dB (nominal Semtech step), with +3 dB as a sensitivity check. The datasheet value of -7.5 dB is used as a second reference.
- **The results are threshold exposure, not a forecast.** They are the share of today's receptions whose best copy lies below an assumed threshold, in a frozen hard-threshold model. A real switch also changes interference, retries and forwarding. Packets that SF8 never decodes are not in the data, so nothing can be said about what SF9 would gain.
- **CoreScope only sees the last hop into an observer.** Repeater-to-repeater margins further up the path are not covered. Observer siting matters a lot, so check the coverage figure.
- **Prefix resolution**: only prefixes that are unique in the catalogue are used for the primary numbers. The share of unresolved copies and of copies without a position is reported.
- **The ALOHA figure is a relative illustration** with assumed load. It ignores LBT, capture, hidden nodes and retries.

## License

MIT, see [LICENSE](LICENSE).
