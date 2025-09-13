# Dpl2LtRtEncoder — Dolby Pro Logic II–style Lt/Rt Encoder (C#)

> **Status:** Hobbyist/educational encoder that approximates Dolby Pro Logic II (PLII) Lt/Rt matrixing. **Not** Dolby‑certified; do not use for broadcast/distribution requiring official compliance.

---

## Overview

**Dpl2LtRtEncoder** converts a 5.1 WAV (WAVEFORMATEXTENSIBLE) into a **stereo Lt/Rt** WAV designed to decode well on **Dolby Pro Logic II** decoders. The implementation follows the classic PLII approach:

* **Center** is mixed into Lt/Rt at \~−3 dB.
* **Surrounds** are matrixed with **asymmetric gains** (near vs. far) and **±90° phase shifts** (via a Hilbert transform):

    * Lt gets **+j** times a weighted sum of the left and right surround feeds.
    * Rt gets **−j** times the *mirrored* weighted sum.
* Surrounds are optionally **high‑pass filtered** and **delayed** (Haas/precedence aid)
* Output is written via a **stereo‑linked look‑ahead limiter** for single‑pass, peak‑safe rendering.

The result is an Lt/Rt that steers dialog to **Center**, keeps music **wide**, and places surround cues to the **rears** on a PLII decoder—close to what you get from a certified encoder for many program types.

---

## Features

* **Input:** 6‑channel WAV (5.1), PCM 16/24/32 or 32‑bit float (Extensible supported)
* **Output:** 2‑channel 32‑bit float WAV (Lt/Rt)
* **Layouts:** `--layout back` (`FL FR FC LFE BL BR`) or `--layout side` (`FL FR FC LFE SL SR`)
* **PLII‑style matrix:** center @ −3 dB (tunable); surrounds with near/far gains (defaults −1.2 dB / −6.2 dB) and ±90° phase via Hilbert FIR
* **Surround pre‑processing:** 2nd‑order Butterworth **high‑pass** (default 120 Hz) and **delay** (default 10 ms)
* **Limiter:** Stereo‑linked **look‑ahead** peak limiter (default 2 ms look‑ahead) with soft knee

---

## Build

```bash
# .NET 8+ recommended
dotnet new console -n Dpl2LtRtEncoder
cd Dpl2LtRtEncoder
dotnet add package NAudio --version 2.2.1
# Replace Program.cs with the one in this repo (the canvas code)
dotnet build
```

---

## Usage

```bash
dotnet run -- <input5_1.wav> <output_ltrt.wav> \
  [--layout back|side] \
  [--sr-delay-ms 10] [--hp-sur-freq 120] \
  [--center-db -3] [--near-sur-db -1.2] [--far-sur-db -6.2] \
  [--lfe-db -inf] [--limit 0.98]
```

### Examples

**Standard 5.1(back) Lt/Rt (BL/BR):**

```bash
dotnet run -- "5.1_in.wav" "out_ltrt.wav" --layout back --center-db -3 --near-sur-db -1.2 --far-sur-db -6.2 --sr-delay-ms 10 --hp-sur-freq 120 --limit 0.98
```


**If your file is 5.1(side):**

```bash
dotnet run -- "5.1_in.wav" "out_ltrt.wav" --layout side
```

Verify channel layout with:

```bash
ffprobe -hide_banner -select_streams a:0 -show_entries stream=channel_layout -of default=nw=1 "5.1_in.wav"
```

---

## Signal Flow

```
   5.1 PCM in  ─┬──────── FL ───────────────┐
                ├──────── FR ───────────────┤
                ├──────── FC ──(×c)─────────┤
                ├──────── LFE ─(×lfe)───────┤      Hilbert(+90°)      Hilbert(−90°)
                ├──────── S1 ─ HPF ─ delay ─┼─→ mix (near S1 + far S2) ──► +j ─┐
                └──────── S2 ─ HPF ─ delay ─┼─→ mix (near S2 + far S1) ──► −j ─┘
                                          Lt sum                      Rt sum
                                               │                          │
                                  limiter (stereo‑linked look‑ahead)      │
                                               │                          │
                                           Lt / Rt (stereo WAV out)
```

Where S1/S2 are **BL/BR** for `--layout back` or **SL/SR** for `--layout side`.

---

## Mathematics

### Channel Matrix (PLII‑style Lt/Rt)

Let inputs be $L, R, C, S_L, S_R$ (LFE optional). Gains:
* $c = 10^{\frac{\text{center-db}}{20}}$  (default −3 dB ⇒ 0.707)
* $a = 10^{\frac{\text{near-sur-db}}{20}}$ (default −1.2 dB ⇒ ≈0.871)
* $b = 10^{\frac{\text{far-sur-db}}{20}}$  (default −6.2 dB ⇒ ≈0.490)
* $\mathcal{H}\{·\}$: Hilbert transform giving approximately **+90°** phase shift (we use FIR windowed ideal)

Then the **Lt/Rt** are:

$$
\begin{aligned}
\text{Lt} &= L + c\,C + j\,\mathcal{H}\{ a\,S_L + b\,S_R \} + \ell \cdot \text{LFE} \\
\text{Rt} &= R + c\,C - j\,\mathcal{H}\{ a\,S_R + b\,S_L \} + \ell \cdot \text{LFE}
\end{aligned}
$$

with optional $\ell$ (LFE gain, default 0). The **asymmetric surround gains** ($a>b$) encourage the decoder’s steering to resolve side/rear direction cleanly and reduce front leakage.

### Surround Pre‑processing

* **High‑pass**: 2nd‑order Butterworth @ `--hp-sur-freq` (default 120 Hz). Implemented with RBJ cookbook biquad coefficients. Purpose: remove LF energy not useful for matrix steering.
* **Delay**: `--sr-delay-ms` (default 10 ms). Psychoacoustic boost for spatial separation (Haas/precedence). Small values (7–15 ms) work well and are conservative for music.

### Hilbert Transform (+90° / −90°)

We approximate a quadrature phase shift using an **odd‑taps FIR Hilbert transformer** (default 201 taps) with a **Blackman window**:

* Ideal impulse response: $h[n] = \frac{2}{\pi n}$ for odd $n\ne 0$, $0$ for even $n$.
* Windowing controls ripple and limits bandwidth; a few hundred taps gives good broadband behavior without audible coloration.
* We compute Lt surround as **+Hilbert** and Rt surround as the **negative** of **+Hilbert** to realize ±90°.

### Look‑ahead Limiter (Stereo‑linked)

We avoid two‑pass normalization in favor of a **single‑pass** limiter that preserves PLII cues:

1. **Look‑ahead buffer** (N samples): store Lt/Rt and track their **max magnitude** over the window.
2. Compute desired gain per frame:

    * $m = \max\{|L|, |R|\}$ over the look‑ahead window
    * $g_d = \min\left(1, \frac{\text{ceiling}}{m + \epsilon}\right)$
3. **Smooth** gain with fast attack, slow release:

    * $g[n] = g_d + (g[n-1] - g_d)\,\alpha$ for attack, or with $\beta$ for release
    * $\alpha = e^{-1/(f_s\,T_{att})},\; \beta = e^{-1/(f_s\,T_{rel})}$
4. Apply **the same** $g[n]$ to Lt and Rt (stereo‑linking) to maintain L/R relationships.

Defaults: look‑ahead 2 ms, attack 0.5 ms, release 80 ms, soft knee, ceiling = `--limit` (e.g., 0.98). For true‑peak safety you can reduce ceiling (e.g., 0.95) or add oversampling (future work).

---

## Parameters & Defaults

* `--layout back|side` — choose 5.1(back) or 5.1(side). Default: `back`.
* `--center-db` — center gain into Lt/Rt (default **−3 dB**).
* `--near-sur-db` — near surround gain (default **−1.2 dB**, ≈0.871).
* `--far-sur-db` — far surround gain (default **−6.2 dB**, ≈0.490).
* `--hp-sur-freq` — surround high‑pass cutoff (default **120 Hz**).
* `--sr-delay-ms` — surround delay in milliseconds (default **10 ms**).
* `--lfe-db` — optional LFE bleed into Lt/Rt (default **−∞** / disabled). If you enable, try **−9 dB**.
* `--limit` — output limiter ceiling (default **0.98**).

---

## Testing & Validation

1. **Channel‑id test:**

    * Play BL‑only and BR‑only clips; verify PLII decoder places them in **rear left** and **rear right**, with minimal center/front leakage.
2. **Dialog anchoring:**

    * Center‑heavy clip should decode with dialog locked to **Center**.
3. **Music imaging:**

    * Wide stereo should remain wide; mono should be stable and centered.
4. **Transient safety:**

    * Check for clipping on loud SFX; limiter should transparently prevent peaks ≥ ceiling.


---

## Known Limitations / Future Work

* **Not Dolby‑certified.** Steering/metadata nuances of official encoders are proprietary.
* **True‑peak limiting** is not implemented; current limiter is sample‑peak based. Add oversampling for TP.
* **Adaptive steering** (program‑dependent behavior) is not modeled; this is a fixed‑matrix approach.
* **Bandwidth of Hilbert** is limited by taps; extremely low or high frequencies may deviate from ±90°.
* **No AC‑3/AAC writer** included. (Potential enhancement: AC‑3 output with *Dolby Surround* flag.)

---

## References & Further Reading

*(Informational, non‑normative; use official Dolby docs for compliance-critical work.)*

* Dolby Surround / Dolby Pro Logic & Pro Logic II overviews (historical background)
* Wikipedia: *Dolby Pro Logic* and *Dolby Pro Logic II* (matrix concepts, ±90° phase for surrounds)
* Robert Bristow‑Johnson, “Audio EQ Cookbook” (biquad formulas)
* Richard G. Lyons, *Understanding Digital Signal Processing* (Hilbert transformers)
* Julius O. Smith III, *Introduction to Digital Filters* (FIR design, windowing)
* Psychoacoustics texts on **precedence (Haas) effect** for spatial localization
* Minnetonka Audio / SurCode product pages (reference for certified workflows)

---

## License & Notice

This encoder is for educational/testing purposes. Dolby, Pro Logic, and related marks are property of Dolby Laboratories. Use at your own risk.
