#!/usr/bin/env python3
"""LoRa time-on-air (Semtech SX126x formula) and theoretical sensitivity estimate.

Usage: lora_airtime.py <payload_bytes> <sf> <bw_hz> <cr 1..4 for 4/5..4/8> [preamble]

Default preamble = MeshCore rule (v1.18.0 RadioLibWrappers.h preambleLengthForSF): 32 symbols for SF <= 8, else 16.
Meshtastic uses 16: pass it explicitly.
"""
import math
import sys

# Demodulator SNR limits per spreading factor (Semtech SX1261/2 datasheet).
SNR_LIMIT_DB = {5: -2.5, 6: -5, 7: -7.5, 8: -10, 9: -12.5, 10: -15, 11: -17.5, 12: -20}


def meshcore_preamble(sf):
    """Preamble length MeshCore sets for a spreading factor (32 for SF <= 8, else 16)."""
    return 32 if sf <= 8 else 16


def airtime_s(payload, sf, bw, cr, preamble=None, crc=True, implicit_header=False):
    if preamble is None:
        preamble = meshcore_preamble(sf)
    # Semtech sx126x_driver sx126x_get_lora_time_on_air_numerator(): SF5/SF6 use 6.25 fixed overhead
    # symbols and no +8 in the ceil numerator; SF>=7 use 4.25 and +8 (LDRO above 16 ms symbol time).
    ts = 2 ** sf / bw
    de = 1 if ts > 0.016 else 0
    num = 8 * payload + 16 * int(crc) - 4 * sf + (0 if implicit_header else 20)
    if sf <= 6:
        den, overhead = 4 * sf, 6.25
    else:
        num += 8
        den, overhead = 4 * (sf - 2 * de), 4.25
    n_payload = 8 + math.ceil(max(num, 0) / den) * (cr + 4)
    return (preamble + overhead) * ts + n_payload * ts


def sensitivity_dbm(sf, bw, noise_figure_db=6.0):
    return -174 + 10 * math.log10(bw) + noise_figure_db + SNR_LIMIT_DB[sf]


if __name__ == "__main__":
    pl, sf, bw, cr = (int(x) for x in sys.argv[1:5])
    pre = int(sys.argv[5]) if len(sys.argv) > 5 else None
    print(f"airtime {airtime_s(pl, sf, bw, cr, pre) * 1000:.1f} ms, sensitivity ~{sensitivity_dbm(sf, bw):.1f} dBm")
