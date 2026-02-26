# Design: elastic_plugin.py Bug Fix

**Date:** 2026-02-26
**Author:** AutoTutorial 3.0
**Scope:** `tools/test_plugins/elastic_plugin.py` only
**Context:** Discovered during Si elastic constants tutorial test (test_report.md in _workspace)

---

## Background

During testing of the elastic constants tutorial, two critical bugs were discovered in `elastic_plugin.py`:

1. `_extract_stru()` regex fails when the tutorial uses `abacustest model inputs` to generate STRU (no inline ATOMIC_SPECIES block exists)
2. `submit_jobs()` only submits a single SCF job — it never runs `abacustest model elastic prepare` to generate the 25 deformed structures required for elastic constant calculation

Root cause of both: the plugin was designed for tutorials that provide STRU inline, but the elastic constants tutorial workflow is:
```
abacustest model inputs → STRU  →  abacustest model elastic prepare → 25 deformed dirs
```

---

## Scope

Only three methods in `elastic_plugin.py` require changes:

| Method | Change |
|--------|--------|
| `prepare_inputs()` | Add STRU fallback via `abacustest model inputs` |
| `submit_jobs()` | Implement full 4-stage pipeline |
| `validate_results()` | Add K-point transparency note |

No changes to `base_plugin.py`, `test_framework_integrated.py`, `relax_plugin.py`, `band_plugin.py`, or `dos_plugin.py`.

---

## Design

### Section 1: `prepare_inputs()` — STRU Fallback

**Problem:** When `_extract_stru()` returns `None` (no inline STRU found), `prepare_inputs()` silently skips STRU creation. The submitted job then lacks STRU, causing ABACUS to crash with `stoi`.

**Solution:** If STRU is missing AND `abacustest_commands` contains `abacustest model inputs`, run it locally to generate the full input directory, then copy the generated files into `relax_dir`.

```
prepare_inputs():
  if stru_content is None:
    if "abacustest model inputs" in abacustest_commands:
      run: abacustest model inputs  (in relax_dir)
      # This generates INPUT, STRU, pseudopotentials, orbitals
    else:
      print("[ERROR] No STRU found and no abacustest model inputs command")
      return []
```

**Fallback priority:**
1. Inline STRU from tutorial → use directly
2. `abacustest model inputs` command found → run to generate
3. Neither → return `[]` with error message

---

### Section 2: `submit_jobs()` — 4-Stage Pipeline

**Problem:** Current `submit_jobs()` submits a single SCF job labeled as "elastic_relax". It never calls `abacustest model elastic prepare` to create 25 deformed structures.

**Solution:** Replace with 4-stage serial pipeline:

```
Stage 1 (existing): submit relax job → wait → download STRU_ION_D

Stage 2 (NEW):      run locally: abacustest model elastic prepare -j elastic_dir
                    # generates org/ + deformed_00 ~ deformed_23 (25 dirs)
                    # copies INPUT, pseudopotentials, orbitals into each

Stage 3 (NEW):      submit all 25 dirs as parallel Bohrium jobs
                    → show dynamic time estimate
                    → ask user to trigger manual check (no auto-polling)
                    → wait for user confirmation → download all 25 results

Stage 4 (NEW):      run locally: abacustest model elastic post -j elastic_dir
                    → parse output → write metrics_elastic.json
```

**Dynamic time estimation (Stage 3):**
Instead of hardcoding "~12 minutes", estimate from job properties:
- n_atoms (from STRU)
- K-points (infer from kspacing or explicit grid)
- n_jobs (25)
- machine_type (from job_config)

Example: Si (8 atoms, 5×5×5 K-points, 25 jobs) → ~10-15 min per job on c32_m64_cpu.

**User check flow (Stage 3):**
```
[预计] 25个任务，Si 8原子，5×5×5 K点，预计每个任务约10-15分钟
[提示] 请等待后运行检查命令：
       python tools/test_framework_integrated.py continue <test_dir>
```
No automatic 30-second polling loop. User decides when to check.

**New helper methods needed:**
- `_prepare_elastic_with_abacustest(elastic_dir)` — runs `abacustest model elastic prepare`
- `_submit_all_deformed_jobs(elastic_dir)` — submits org + deformed_00~23
- `_wait_for_all_jobs(job_ids)` — shows estimate + waits for user signal
- `_post_process_elastic(elastic_dir)` — runs `abacustest model elastic post`, writes metrics_elastic.json
- `_estimate_job_time(stru_path, kspacing, n_jobs, machine_type)` — dynamic estimation

---

### Section 3: `validate_results()` — K-point Transparency

**Problem:** When actual K-points differ from tutorial's expected K-points, results show >5% error with no explanation. The user sees "FAIL" without understanding why.

**Solution:** When loading `metrics_elastic.json`, check for kspacing/K-grid discrepancy. If found, add a note to `validation.warnings`:

```
⚠️ K点说明：实际使用 kspacing=0.14 → [5,5,5]，教程预期值基于 8×8×8。
   数值差异属于正常的K点收敛效应，不影响方法论正确性。
```

Also update `generate_report_section()` to display this warning prominently when comparisons partially fail.

---

## Implementation Plan

1. Add `_prepare_elastic_with_abacustest()` method
2. Add `_submit_all_deformed_jobs()` method
3. Add `_estimate_job_time()` method
4. Add `_wait_for_all_jobs()` method (show estimate + user-triggered check)
5. Add `_post_process_elastic()` method
6. Modify `prepare_inputs()` — add STRU fallback branch
7. Modify `submit_jobs()` — replace single SCF with 4-stage pipeline
8. Modify `validate_results()` — add K-point discrepancy warning

---

## Test Verification

After implementation, run against:
- `_workspace/20260203_105918_弹性常数计算 copy 3/07_final.md` (Si + TiO₂)

Expected behavior:
1. `prepare_inputs()`: STRU fallback triggers, generates full input via abacustest
2. `submit_jobs()`: All 4 stages execute, 25 jobs submitted, metrics_elastic.json created
3. `validate_results()`: K-point note appears when kspacing=0.14 used

---

## Non-Goals

- No changes to `base_plugin.py` TestInfo/ValidationResult dataclasses
- No changes to `relax_plugin.py`, `band_plugin.py`, `dos_plugin.py`
- No changes to the 4-stage pattern itself (it's already in submit_jobs conceptually)
- No changes to job submission infrastructure (BohriumJobManager)
