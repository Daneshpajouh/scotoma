# Checking the C_9 to C_12 duals

These four files are the rational duals behind the floors C_9 = 9, C_10 = 7, C_11 = 6 and C_12 = 4 at spacer length 20 (word length N = 19). Each is a nonnegative measure on the target set. If every radius-m ball carries load at most 1 and the total mass is above r - 1, every cover needs at least r balls.

## Inputs (SHA-256)

    0f3ce51c8fd6cbd406c6866de1a5c8b948cbf0dacc952277d4992f9162118a95  data/floor_certificates/cert_m9_r9.json
    4577123ca0157bcbfbb7bf9cf3f46301488d4123983b2f02e7f6d7d4c1426370  data/floor_certificates/cert_m10_r7.json
    9d8ae33019f45ef9cb6a0b4149c29ecae438e2d7d0195f1b294b3f17cd4d12b2  data/floor_certificates/cert_m11_r6.json
    896515839a68a2358d17f1385efe6fef116c3de4ce265b0b79240494c34b37b9  data/floor_certificates/cert_m12_r4.json
    8ae5295da593f4d2e15fab680e886ac5e2bf84ab3040d97d2ba0bb6b8bd752da  cert_19_9_9.json
    1b658bbe65ff868ea0a2f32617b6dc8d02f69378bee776930c2bd9410703a91d  cert_19_10_7.json
    cde4688e8af73f30ef8254e8e263971cb687f8562c6347cb3e845019a6a785fc  cert_19_11_6.json

The first four are the released duals. The last three are the solver's unnormalised output for m = 9, 10 and 11.

## 1. The exhaustive check (the one the proof rests on)

Python 3, standard library only, about one minute:

    python3 -I verify_ceiling_exhaustive.py data/floor_certificates/cert_m9_r9.json \
        data/floor_certificates/cert_m10_r7.json data/floor_certificates/cert_m11_r6.json \
        data/floor_certificates/cert_m12_r4.json

Expected: exit code 0, one line per file, then `ALL PASS`. The key values on each line:

    cert_m9_r9.json   mass 8999968/1104275   78732 of 78732 classes swept   max load 1 at prefix 10011000 free weight 0
    cert_m10_r7.json  mass 13999966/2124163  10206 of 10206 classes swept   max load 1 at prefix 1X1110 free weight 0
    cert_m11_r6.json  mass 2999993/582027    3645 of 3645 classes swept     max load 1 at prefix 10011 free weight 0
    cert_m12_r4.json  mass 7843719/2500000   459 of 459 classes swept       max load 9999999/10000000 at prefix 000 free weight 0

Each mass is above r - 1 (8, 6, 5 and 3). Any failure prints a line starting `FAIL:` and exits 1.

What it checks: the format of every field the proof uses (masses must be exact: JSON strings or integers; decimal numbers are decoded exactly, never through a float); nonnegative masses that sum to the stated total; every block within distance m of the library; every centre class at load at most 1, with the maximum equal to the stated one; the named maximising class, if the file names one; on 40 sampled classes per dual, a random centre of the class against the load the sweep recorded. Descriptive fields, such as `centre_classes_verified`, are printed and not checked. At m = 12 the solver swept the 136 binary-prefix classes by design; this check sweeps all 459.

## 2. From solver output to released file

    python3 -I emit_ceiling_certificate.py cert_19_9_9.json

This writes data/floor_certificates/cert_m9_r9.json. Run it in a copy of the release folder, because it overwrites that file and updates INDEX.json. It first runs verify_F_independent.py and verify_ceiling_certificate.py and stops if either fails. Done for all three solver files, the output is byte for byte the released cert_m9_r9.json, cert_m10_r7.json and cert_m11_r6.json (hashes above).

## 3. The sampled consistency check

    python3 -I verify_ceiling_certificate.py cert_19_9_9.json

This checks the arithmetic, equal loads on six sampled pairs, and the load at the library centres and 90 sampled centres. It does not establish the maximum load, so it does not prove the floor on its own. Use check 1 for that.
