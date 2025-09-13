// Dolby Pro Logic II Lt/Rt Encoder (Approximate, non‑certified)
// -------------------------------------------------------------
// Console app that reads a 5.1 WAV (WAVEFORMATEXTENSIBLE) and writes a stereo WAV
// matrix-encoded for playback with Dolby Pro Logic II decoders.
//
// Features:
//  - PLII-style asymmetric surround weighting (near/far)
//  - Surround high-pass (Butterworth, 2nd-order, default 120 Hz)
//  - Optional surround delay (default 10 ms)
//  - Hilbert-transform (FIR, windowed) to apply ±90° phase to surround sums
//  - Optional tiny LFE bleed (default 0)
//  - Peak-safe output normalization to a chosen ceiling (default 0.98)
//  - Channel layout chosen via flag: 5.1(back) or 5.1(side) + optional surround swap
//
// DISCLAIMER: This approximates Pro Logic II encoding and is NOT Dolby-certified.
//             For broadcast/commercial use, use an official encoder (e.g., SurCode).
//
// Build:
//   dotnet new console -n Dpl2LtRtEncoder
//   cd Dpl2LtRtEncoder
//   dotnet add package NAudio --version 2.2.1
//   Replace Program.cs with this file content
//   dotnet build
//
// Usage:
//   dotnet run -- <input.wav> <output.wav> [--layout back|side] [--sr-delay-ms 10] \
//               [--hp-sur-freq 120] [--center-db -3] [--near-sur-db -1.2] [--far-sur-db -6.2] \
//               [--lfe-db -inf] [--limit 0.98]
//
// Example (file tagged 5.1(back) with BL/BR):
//   dotnet run -- 5.1_in.wav out_ltrt.wav --layout back --center-db -3 --near-sur-db -1.2 --far-sur-db -6.2

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using NAudio.Wave;

namespace Dpl2LtRtEncoder
{
    class Program
    {
        static int Main(string[] args)
        {
            if (args.Length < 2)
            {
                Console.WriteLine("Usage: Dpl2LtRtEncoder <input.wav> <output.wav> [--layout back|side] [--sr-delay-ms 10] [--hp-sur-freq 120] [--center-db -3] [--near-sur-db -1.2] [--far-sur-db -6.2] [--lfe-db -inf] [--limit 0.98]");
                return 1;
            }

            var inputPath = args[0];
            var outputPath = args[1];

            var layout = "back"; // default; you can pass --layout side
            float surroundDelayMs = 10f;
            float hpSurFreq = 120f;
            float centerDb = -3f;
            float nearSurDb = -1.2f; // ~0.871
            float farSurDb = -6.2f;  // ~0.490
            float lfeDb = float.NegativeInfinity; // default: no LFE
            float limitCeiling = 0.98f;

            // Parse simple args
            for (int i = 2; i < args.Length; i++)
            {
                switch (args[i])
                {
                    case "--layout": layout = args[++i]; break;
                    case "--sr-delay-ms": surroundDelayMs = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--hp-sur-freq": hpSurFreq = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--center-db": centerDb = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--near-sur-db": nearSurDb = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--far-sur-db": farSurDb = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--lfe-db": lfeDb = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                    case "--limit": limitCeiling = float.Parse(args[++i], CultureInfo.InvariantCulture); break;
                }
            }

            if (!File.Exists(inputPath))
            {
                Console.WriteLine($"Input not found: {inputPath}");
                return 1;
            }

            try
            {
                using var reader = new WaveFileReader(inputPath);
                var wf = reader.WaveFormat;
                if (wf.Channels != 6)
                {
                    Console.WriteLine($"Error: input must be 6 channels (found {wf.Channels}).");
                    return 1;
                }

                int sampleRate = wf.SampleRate;

                // Choose channel mapping (assumes interleaved order matches layout choice)
                ChannelMap mapping = layout.Equals("side", StringComparison.OrdinalIgnoreCase)
                    ? ChannelMapper.FivePointOneSide()
                    : ChannelMapper.FivePointOneBack();

                Console.WriteLine($"SampleRate: {sampleRate} Hz, Bits: {wf.BitsPerSample}, Channels: {wf.Channels}");
                Console.WriteLine($"Using layout: {mapping.LayoutName}");

                // Input format details
                int channels = wf.Channels;
                int bytesPerSample = wf.BitsPerSample / 8;
                int bytesPerFrame = wf.BlockAlign; // bytes per multi-channel frame

                // Create writer and limiter (streaming)
                var outFormat = WaveFormat.CreateIeeeFloatWaveFormat(sampleRate, 2);
                using var writer = new WaveFileWriter(outputPath, outFormat);
                var limiter = new LookaheadLimiterStereo(sampleRate, lookaheadMs: 2.0f, attackMs: 0.5f, releaseMs: 80f, ceiling: limitCeiling, kneeDb: 3.0f);
                float[] frameOut = new float[2];

                // DSP blocks
                var hpBL = new BiquadHighPass(sampleRate, hpSurFreq, 0.707f);
                var hpBR = new BiquadHighPass(sampleRate, hpSurFreq, 0.707f);
                int hilbertTaps = 201; // odd; trade-off between accuracy and CPU/latency
                var hilbertLt = new HilbertTransformer(hilbertTaps); // +90° for Lt surround sum
                var hilbertRt = new HilbertTransformer(hilbertTaps); // +90° for Rt surround sum (we negate later for -90°)
                int delaySamples = (int)Math.Round(surroundDelayMs * sampleRate / 1000.0);
                var delayBL = new DelayLine(delaySamples);
                var delayBR = new DelayLine(delaySamples);

                float cGain = DbToLin(centerDb);
                float nearGain = DbToLin(nearSurDb);
                float farGain = DbToLin(farSurDb);
                float lfeGain = float.IsNegativeInfinity(lfeDb) ? 0f : DbToLin(lfeDb);

                // Read loop: multiple of frames per buffer
                int framesPerBuffer = 1024;
                byte[] buf = new byte[bytesPerFrame * framesPerBuffer];

                while (true)
                {
                    int bytesRead = reader.Read(buf, 0, buf.Length);
                    if (bytesRead <= 0) break;
                    int frames = bytesRead / bytesPerFrame;

                    for (int n = 0; n < frames; n++)
                    {
                        int baseOffset = n * bytesPerFrame;

                        // Extract per-channel samples to float [-1,1]
                        float[] s = new float[channels];
                        for (int ch = 0; ch < channels; ch++)
                        {
                            int off = baseOffset + ch * bytesPerSample;
                            switch (wf.Encoding)
                            {
                                case WaveFormatEncoding.Pcm:
                                    if (wf.BitsPerSample == 16)
                                    {
                                        short v = (short)(buf[off] | (buf[off + 1] << 8));
                                        s[ch] = Math.Clamp(v / 32768f, -1f, 1f);
                                    }
                                    else if (wf.BitsPerSample == 24)
                                    {
                                        int v = buf[off] | (buf[off + 1] << 8) | (buf[off + 2] << 16);
                                        if ((v & 0x800000) != 0) v |= unchecked((int)0xFF000000); // sign extend
                                        s[ch] = Math.Clamp(v / 8388608f, -1f, 1f);
                                    }
                                    else if (wf.BitsPerSample == 32)
                                    {
                                        int v = BitConverter.ToInt32(buf, off);
                                        s[ch] = Math.Clamp(v / 2147483648f, -1f, 1f);
                                    }
                                    else
                                    {
                                        throw new InvalidOperationException($"Unsupported PCM bit depth: {wf.BitsPerSample}");
                                    }
                                    break;
                                case WaveFormatEncoding.IeeeFloat:
                                    if (wf.BitsPerSample == 32)
                                    {
                                        s[ch] = BitConverter.ToSingle(buf, off);
                                    }
                                    else
                                    {
                                        throw new InvalidOperationException($"Unsupported IEEE float bit depth: {wf.BitsPerSample}");
                                    }
                                    break;
                                case WaveFormatEncoding.Extensible:
                                    // WAVEFORMATEXTENSIBLE often still contains PCM/Float data; infer by BitsPerSample
                                    if (wf.BitsPerSample == 16)
                                    {
                                        short v16 = (short)(buf[off] | (buf[off + 1] << 8));
                                        s[ch] = Math.Clamp(v16 / 32768f, -1f, 1f);
                                    }
                                    else if (wf.BitsPerSample == 24)
                                    {
                                        int v = buf[off] | (buf[off + 1] << 8) | (buf[off + 2] << 16);
                                        if ((v & 0x800000) != 0) v |= unchecked((int)0xFF000000);
                                        s[ch] = Math.Clamp(v / 8388608f, -1f, 1f);
                                    }
                                    else if (wf.BitsPerSample == 32)
                                    {
                                        // Heuristic: assume float if valid float range, else PCM32
                                        float f = BitConverter.ToSingle(buf, off);
                                        if (float.IsFinite(f) && Math.Abs(f) <= 4f)
                                            s[ch] = Math.Clamp(f, -1f, 1f);
                                        else
                                        {
                                            int vi = BitConverter.ToInt32(buf, off);
                                            s[ch] = Math.Clamp(vi / 2147483648f, -1f, 1f);
                                        }
                                    }
                                    else
                                    {
                                        throw new InvalidOperationException($"Unsupported Extensible bit depth: {wf.BitsPerSample}");
                                    }
                                    break;
                                default:
                                    throw new InvalidOperationException($"Unsupported source encoding: {wf.Encoding}");
                            }
                        }

                        // Map channels (indices assume interleaved order according to chosen layout)
                        float FL = s[mapping.FL];
                        float FR = s[mapping.FR];
                        float FC = s[mapping.FC];
                        float LFE = s[mapping.LFE];
                        float S1 = s[mapping.SL_or_BL]; // SL or BL
                        float S2 = s[mapping.SR_or_BR]; // SR or BR

                        // Surround preprocessing: HPF -> delay
                        float S1p = delayBL.Process(hpBL.Process(S1));
                        float S2p = delayBR.Process(hpBR.Process(S2));

                        // PLII-style surround contributions (no pre-negation)
                        // Lt: +j*(near*BL/SL + far*BR/SR)
                        // Rt: -j*(near*BR/SR + far*BL/SL)
                        float S_lt_raw = nearGain * S1p + farGain * S2p;   // Lt near=BL/SL, far=BR/SR
                        float S_rt_raw = nearGain * S2p + farGain * S1p;   // Rt near=BR/SR, far=BL/SL

                        // Hilbert ±90°
                        float S_lt = hilbertLt.Process(S_lt_raw);   // +90°
                        float S_rt = -hilbertRt.Process(S_rt_raw);  // -90°

                        float Lt = FL + cGain * FC + S_lt + lfeGain * LFE;
                        float Rt = FR + cGain * FC + S_rt + lfeGain * LFE;

                        if (limiter.Process(Lt, Rt, out float yL, out float yR))
                        {
                            frameOut[0] = yL; frameOut[1] = yR;
                            writer.WriteSamples(frameOut, 0, 2);
                        }
                    }
                }

                // Flush limiter tail
                while (limiter.Flush(out float yL2, out float yR2))
                {
                    frameOut[0] = yL2; frameOut[1] = yR2;
                    writer.WriteSamples(frameOut, 0, 2);
                }

                Console.WriteLine($"Done. Wrote Lt/Rt to: {outputPath}");
                return 0;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error: " + ex.Message);
                return 1;
            }
        }

        static float DbToLin(float db) => (float)Math.Pow(10.0, db / 20.0);

        static void NormalizeTwoPass(List<float> L, List<float> R, float ceiling)
        {
            float peak = 1e-9f;
            int count = L.Count;
            for (int i = 0; i < count; i++)
            {
                float m = Math.Abs(L[i]);
                if (m > peak) peak = m;
                m = Math.Abs(R[i]);
                if (m > peak) peak = m;
            }
            if (peak <= ceiling) return; // already under
            float g = ceiling / peak;
            for (int i = 0; i < count; i++)
            {
                L[i] *= g;
                R[i] *= g;
            }
        }
    }

    // Maps channel indices for common 5.1 layouts
    class ChannelMap
    {
        public int FL, FR, FC, LFE, SL_or_BL, SR_or_BR;
        public string LayoutName;
        public ChannelMap(int fl, int fr, int fc, int lfe, int s1, int s2, string name)
        { FL = fl; FR = fr; FC = fc; LFE = lfe; SL_or_BL = s1; SR_or_BR = s2; LayoutName = name; }
    }

    static class ChannelMapper
    {
        // 5.1(back): FL FR FC LFE BL BR (assumed interleaved order)
        public static ChannelMap FivePointOneBack() => new ChannelMap(0, 1, 2, 3, 4, 5, "5.1(back): FL FR FC LFE BL BR");
        // 5.1(side): FL FR FC LFE SL SR (assumed interleaved order)
        public static ChannelMap FivePointOneSide() => new ChannelMap(0, 1, 2, 3, 4, 5, "5.1(side): FL FR FC LFE SL SR");
    }

    // Simple delay line (samples)
    class DelayLine
    {
        private readonly float[] _buf;
        private int _w, _r;
        public DelayLine(int samples)
        {
            if (samples < 0) samples = 0;
            _buf = new float[Math.Max(1, samples)];
            _w = 0; _r = 0;
        }
        public float Process(float x)
        {
            float y = _buf[_r];
            _buf[_w] = x;
            _w++; if (_w >= _buf.Length) _w = 0;
            _r++; if (_r >= _buf.Length) _r = 0;
            return y;
        }
    }

    // 2nd-order Butterworth High-Pass (Biquad)
    class BiquadHighPass
    {
        private readonly float fb0, fb1, fb2; // feedforward (normalized)
        private readonly float fa1, fa2;      // feedback (normalized, a0=1)
        private float z1, z2;
        public BiquadHighPass(int fs, float fc, float q)
        {
            if (fc < 10) fc = 10; // avoid DC weirdness
            float w0 = 2.0f * (float)Math.PI * (fc / fs);
            float cosw0 = (float)Math.Cos(w0);
            float sinw0 = (float)Math.Sin(w0);
            float alpha = sinw0 / (2.0f * q);

            // RBJ cookbook high-pass
            float b0 = (1 + cosw0) / 2.0f;
            float b1 = -(1 + cosw0);
            float b2 = (1 + cosw0) / 2.0f;
            float a0 = 1 + alpha;
            float a1 = -2 * cosw0;
            float a2 = 1 - alpha;

            fb0 = b0 / a0;
            fb1 = b1 / a0;
            fb2 = b2 / a0;
            fa1 = a1 / a0;
            fa2 = a2 / a0;
        }
        public float Process(float x)
        {
            // Transposed Direct Form II
            float y = fb0 * x + z1;
            z1 = fb1 * x - fa1 * y + z2;
            z2 = fb2 * x - fa2 * y;
            return y;
        }
    }

    // FIR Hilbert transformer (+90°) using windowed ideal impulse response
    class HilbertTransformer
    {
        private readonly float[] h; // taps
        private readonly float[] d; // delay line
        private int idx;

        public HilbertTransformer(int taps)
        {
            if (taps % 2 == 0) taps += 1; // ensure odd
            h = MakeHilbertTaps(taps);
            d = new float[taps];
            idx = 0;
        }

        public float Process(float x)
        {
            d[idx] = x;
            float acc = 0f;
            int j = idx;
            for (int i = 0; i < h.Length; i++)
            {
                acc += h[i] * d[j];
                if (--j < 0) j = d.Length - 1;
            }

            if (++idx >= d.Length) idx = 0;
            return acc;
        }

        // Ideal Hilbert impulse response: h[n] = 2/(pi*n) for n odd, 0 for n even, n!=0
        // Windowed with Blackman to control ripple.
        private static float[] MakeHilbertTaps(int taps)
        {
            int M = taps;
            int mid = M / 2;
            var coeffs = new float[M];
            for (int n = 0; n < M; n++)
            {
                int k = n - mid;
                float hn = 0f;
                if (k != 0 && (k % 2 != 0)) // odd k
                {
                    hn = (float)(2.0 / (Math.PI * k)); // ideal
                }

                // Blackman window
                float w = 0.42f - 0.5f * (float)Math.Cos(2 * Math.PI * n / (M - 1)) +
                          0.08f * (float)Math.Cos(4 * Math.PI * n / (M - 1));
                coeffs[n] = hn * w;
            }

            return coeffs;
        }
        // Stereo-linked look-ahead peak limiter (sample-peak, with headroom). For true-peak safety, set ceiling ≤ 0.98.

    }

    class LookaheadLimiterStereo
    {
        private readonly int N;
        private readonly float ceiling;
        private readonly float attCoeff, relCoeff;
        private readonly float[] bufL, bufR, magBuf;
        private int w = 0, r = 0, count = 0;
        private float env = 1f;

        public LookaheadLimiterStereo(int sampleRate, float lookaheadMs = 2.0f, float attackMs = 0.5f, float releaseMs = 80f, float ceiling = 0.98f, float kneeDb = 3.0f)
        {
            N = Math.Max(1, (int)Math.Round(sampleRate * lookaheadMs / 1000.0));
            this.ceiling = ceiling;
            attCoeff = (float)Math.Exp(-1.0 / (Math.Max(1.0, attackMs) * sampleRate / 1000.0));
            relCoeff = (float)Math.Exp(-1.0 / (Math.Max(1.0, releaseMs) * sampleRate / 1000.0));
            bufL = new float[N];
            bufR = new float[N];
            magBuf = new float[N];
        }

        public bool Process(float inL, float inR, out float outL, out float outR)
        {
            // Write incoming sample
            float m = MathF.Max(MathF.Abs(inL), MathF.Abs(inR));
            bufL[w] = inL; bufR[w] = inR; magBuf[w] = m;
            w++; if (w >= N) w = 0;
            if (count < N) { count++; outL = 0; outR = 0; return false; }

            // Output oldest sample at r
            float xL = bufL[r];
            float xR = bufR[r];

            // Compute max over look-ahead window (simple scan, N is small ~96 @ 48kHz 2ms)
            float aheadMax = 0f;
            for (int i = 0; i < N; i++)
            {
                float mm = magBuf[i];
                if (mm > aheadMax) aheadMax = mm;
            }

            float desired = 1f;
            if (aheadMax > ceiling && aheadMax > 1e-12f)
                desired = ceiling / aheadMax;

            // Smooth gain (attack fast, release slow)
            if (desired < env)
                env = desired + (env - desired) * attCoeff;
            else
                env = desired + (env - desired) * relCoeff;

            outL = xL * env;
            outR = xR * env;

            // Advance read index
            r++; if (r >= N) r = 0;
            return true;
        }

        public bool Flush(out float outL, out float outR)
        {
            if (count == 0) { outL = 0; outR = 0; return false; }

            // Compute max over remaining window
            float aheadMax = 0f;
            for (int i = 0; i < count; i++)
            {
                int idx = (r + i) % N;
                float mm = magBuf[idx];
                if (mm > aheadMax) aheadMax = mm;
            }
            float desired = 1f;
            if (aheadMax > ceiling && aheadMax > 1e-12f)
                desired = ceiling / aheadMax;
            if (desired < env)
                env = desired + (env - desired) * attCoeff;
            else
                env = desired + (env - desired) * relCoeff;

            outL = bufL[r] * env;
            outR = bufR[r] * env;

            // Clear this sample and advance
            magBuf[r] = 0; bufL[r] = 0; bufR[r] = 0;
            r++; if (r >= N) r = 0;
            count--;
            return true;
        }
    }
}
