# Python Bug Fix Plan: Duplicate Lines and Lost Candidates

## Overview

This document provides a systematic plan to fix two bugs in the Python `rgmatch` implementation (`src/rgmatch/core.py`):

1. **Bug #1**: 99 duplicate lines in output (same line appearing twice)
2. **Bug #2**: 62 lost proximity candidates (biologically relevant associations silently discarded)

Both bugs stem from the same root cause: misuse of `exon_down` and `exon_up` variables.

---

## Root Cause

### Current (Buggy) Behavior

The variables `exon_down` and `exon_up` serve **two conflicting purposes**:

1. **Purpose A**: Track the closest non-overlapping (proximity) candidate for final output
2. **Purpose B**: Temporarily hold overlapping candidates before immediate output

This dual-use causes:
- **Duplicates**: Overlapping candidates are output immediately AND again at the end
- **Lost candidates**: Proximity candidates are overwritten by overlapping candidates

### Code Location

**File**: `src/rgmatch/core.py`

**Problematic sections**:
- Lines 1093-1111 (Case 2: Partial overlap left)
- Lines 1153-1173 (Case 3: Exon inside region - upstream portion)
- Lines 1191-1210 (Case 3: Exon inside region - downstream portion)
- Lines 1250-1270 (Case 4: Partial overlap right)
- Lines 1335-1347 (Final output section)

---

## The Fix

### Strategy

**Separate the two purposes** by NOT using `exon_down`/`exon_up` for overlapping candidates. Instead, output overlapping candidates directly without storing them in these variables.

### Detailed Changes

#### Change 1: Case 2 - Partial Overlap Left (lines 1093-1111)

**Current code (buggy)**:
```python
if exon.getEnd() < end:
    if isLastExon is True:
        region_overlap = end - exon.getEnd()
        pctg_region = (float(region_overlap) / region_length) * 100
        if mygene.getStrand() == "+":
            tag = "DOWNSTREAM"
            exon_down = Candidate(...)  # BUG: Overwrites proximity candidate!
            if tts > 0:
                mychecks = checkTTS(start, end, exon_down)
                for assoc in mychecks:
                    myfinaloutput.append(Candidate(...))
            else:
                myfinaloutput.append(exon_down)
        else:
            tag = "UPSTREAM"
            exon_up = Candidate(...)  # BUG: Overwrites proximity candidate!
            mychecks = checkTSS(start, end, exon_up)
            for assoc in mychecks:
                myfinaloutput.append(Candidate(...))
```

**Fixed code**:
```python
if exon.getEnd() < end:
    if isLastExon is True:
        region_overlap = end - exon.getEnd()
        pctg_region = (float(region_overlap) / region_length) * 100
        if mygene.getStrand() == "+":
            tag = "DOWNSTREAM"
            # FIX: Use local variable, don't overwrite exon_down
            overlap_candidate = Candidate(exon.getStart(), exon.getEnd(), mygene.getStrand(),
                exon.getExon(), tag, mytranscript.getTranscriptID(), mygene.getGeneID(),
                0, pctg_region, -1, TSSdistance)
            if tts > 0:
                mychecks = checkTTS(start, end, overlap_candidate)
                for assoc in mychecks:
                    myfinaloutput.append(Candidate(overlap_candidate.getStart(),
                        overlap_candidate.getEnd(), overlap_candidate.getStrand(),
                        overlap_candidate.getExonNr(), assoc[0], overlap_candidate.getTranscript(),
                        overlap_candidate.getGene(), overlap_candidate.getDistance(),
                        assoc[1], assoc[2], TSSdistance))
            else:
                myfinaloutput.append(overlap_candidate)
        else:
            tag = "UPSTREAM"
            # FIX: Use local variable, don't overwrite exon_up
            overlap_candidate = Candidate(exon.getStart(), exon.getEnd(), mygene.getStrand(),
                exon.getExon(), tag, mytranscript.getTranscriptID(), mygene.getGeneID(),
                0, pctg_region, -1, TSSdistance)
            mychecks = checkTSS(start, end, overlap_candidate)
            for assoc in mychecks:
                myfinaloutput.append(Candidate(overlap_candidate.getStart(),
                    overlap_candidate.getEnd(), overlap_candidate.getStrand(),
                    overlap_candidate.getExonNr(), assoc[0], overlap_candidate.getTranscript(),
                    overlap_candidate.getGene(), overlap_candidate.getDistance(),
                    assoc[1], assoc[2], TSSdistance))
```

#### Change 2: Case 3 - Exon Inside Region, Upstream Portion (lines 1153-1173)

**Current code (buggy)**:
```python
if start < exon.getStart():
    if isFirstExon is True:
        region_overlap = exon.getStart() - start
        pctg_region = (float(region_overlap) / region_length) * 100

        if mygene.getStrand() == "-":
            tag = "DOWNSTREAM"
            exon_down = Candidate(...)  # BUG!
            if tts > 0:
                mychecks = checkTTS(start, end, exon_down)
                for assoc in mychecks:
                    myfinaloutput.append(Candidate(...))
            else:
                myfinaloutput.append(exon_down)

        else:
            tag = "UPSTREAM"
            exon_up = Candidate(...)  # BUG!
            mychecks = checkTSS(start, end, exon_up)
            for assoc in mychecks:
                myfinaloutput.append(Candidate(...))
```

**Fixed code**:
```python
if start < exon.getStart():
    if isFirstExon is True:
        region_overlap = exon.getStart() - start
        pctg_region = (float(region_overlap) / region_length) * 100

        if mygene.getStrand() == "-":
            tag = "DOWNSTREAM"
            # FIX: Use local variable
            overlap_candidate = Candidate(exon.getStart(), exon.getEnd(), mygene.getStrand(),
                exon.getExon(), tag, mytranscript.getTranscriptID(), mygene.getGeneID(),
                0, pctg_region, -1, TSSdistance)
            if tts > 0:
                mychecks = checkTTS(start, end, overlap_candidate)
                for assoc in mychecks:
                    myfinaloutput.append(Candidate(overlap_candidate.getStart(),
                        overlap_candidate.getEnd(), overlap_candidate.getStrand(),
                        overlap_candidate.getExonNr(), assoc[0], overlap_candidate.getTranscript(),
                        overlap_candidate.getGene(), overlap_candidate.getDistance(),
                        assoc[1], assoc[2], TSSdistance))
            else:
                myfinaloutput.append(overlap_candidate)

        else:
            tag = "UPSTREAM"
            # FIX: Use local variable
            overlap_candidate = Candidate(exon.getStart(), exon.getEnd(), mygene.getStrand(),
                exon.getExon(), tag, mytranscript.getTranscriptID(), mygene.getGeneID(),
                0, pctg_region, -1, TSSdistance)
            mychecks = checkTSS(start, end, overlap_candidate)
            for assoc in mychecks:
                myfinaloutput.append(Candidate(overlap_candidate.getStart(),
                    overlap_candidate.getEnd(), overlap_candidate.getStrand(),
                    overlap_candidate.getExonNr(), assoc[0], overlap_candidate.getTranscript(),
                    overlap_candidate.getGene(), overlap_candidate.getDistance(),
                    assoc[1], assoc[2], TSSdistance))
```

#### Change 3: Case 3 - Exon Inside Region, Downstream Portion (lines 1191-1210)

Apply the same pattern: replace `exon_down = Candidate(...)` and `exon_up = Candidate(...)` with `overlap_candidate = Candidate(...)`.

#### Change 4: Case 4 - Partial Overlap Right (lines 1250-1270)

Apply the same pattern: replace `exon_down = Candidate(...)` and `exon_up = Candidate(...)` with `overlap_candidate = Candidate(...)`.

---

## Summary of All Changes

| Location | Line Numbers | Change |
|----------|--------------|--------|
| Case 2 (positive strand) | ~1099 | `exon_down = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 2 (negative strand) | ~1108 | `exon_up = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 3 upstream (negative strand) | ~1160 | `exon_down = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 3 upstream (positive strand) | ~1170 | `exon_up = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 3 downstream (positive strand) | ~1198 | `exon_down = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 3 downstream (negative strand) | ~1207 | `exon_up = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 4 (negative strand) | ~1257 | `exon_down = Candidate(...)` → `overlap_candidate = Candidate(...)` |
| Case 4 (positive strand) | ~1267 | `exon_up = Candidate(...)` → `overlap_candidate = Candidate(...)` |

**Total**: 8 variable assignments need to be changed from `exon_down`/`exon_up` to a local `overlap_candidate`.

---

## Verification

After applying the fix, run:

```bash
# Generate new output
python -m rgmatch -g test_data/benchmark.gtf -b test_data/benchmark.bed -o /tmp/fixed_output.txt

# Check for duplicates (should be 0)
sort /tmp/fixed_output.txt | uniq -d | wc -l
# Expected: 0

# Check total unique lines (should be 15,950 = 15,888 + 62 recovered)
sort /tmp/fixed_output.txt | uniq | wc -l
# Expected: 15,950

# Compare with Rust output (should match exactly)
sort /tmp/fixed_output.txt | uniq > /tmp/fixed_uniq.txt
sort /tmp/rust_output.txt | uniq > /tmp/rust_uniq.txt
diff /tmp/fixed_uniq.txt /tmp/rust_uniq.txt
# Expected: no differences
```

---

## Expected Results After Fix

| Metric | Before Fix | After Fix |
|--------|------------|-----------|
| Total lines | 15,987 | 15,950 |
| Unique lines | 15,888 | 15,950 |
| Duplicate lines | 99 | 0 |
| Lost candidates | 62 | 0 |

---

## Biological Impact

The fix will:

1. **Eliminate false duplicates**: No more inflated output with repeated lines
2. **Recover lost associations**: 62 biologically relevant region-gene associations will be restored
3. **Improve reproducibility**: Results will no longer depend on gene processing order
4. **Match comprehensive Rust output**: Both implementations will produce identical results

---

## Alternative Approach (More Invasive)

If a cleaner refactor is desired, consider separating the tracking entirely:

```python
# Instead of single variables:
exon_down = None
exon_up = None

# Use separate tracking for proximity vs overlap:
proximity_down = None  # For Case 1/6 (non-overlapping)
proximity_up = None    # For Case 1/6 (non-overlapping)
# Overlapping candidates go directly to myfinaloutput
```

This makes the intent clearer but requires more extensive changes.

---

## Files to Modify

1. **Primary**: `src/rgmatch/core.py`
   - Lines 1093-1111 (Case 2)
   - Lines 1153-1173 (Case 3 upstream)
   - Lines 1191-1210 (Case 3 downstream)
   - Lines 1250-1270 (Case 4)

2. **Tests**: Update golden output files to reflect the corrected behavior
   - `test_data/golden_output.txt` (will have 15,950 unique lines instead of 15,888)

---

## Implementation Checklist

- [ ] Change Case 2 positive strand: `exon_down` → `overlap_candidate`
- [ ] Change Case 2 negative strand: `exon_up` → `overlap_candidate`
- [ ] Change Case 3 upstream negative strand: `exon_down` → `overlap_candidate`
- [ ] Change Case 3 upstream positive strand: `exon_up` → `overlap_candidate`
- [ ] Change Case 3 downstream positive strand: `exon_down` → `overlap_candidate`
- [ ] Change Case 3 downstream negative strand: `exon_up` → `overlap_candidate`
- [ ] Change Case 4 negative strand: `exon_down` → `overlap_candidate`
- [ ] Change Case 4 positive strand: `exon_up` → `overlap_candidate`
- [ ] Run verification tests
- [ ] Update golden output files
- [ ] Update Rust implementation to remove the "Python compatibility" fix (optional)
