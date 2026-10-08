# Singleton Extension in the BWA FMD Index

## Overview

`BwaFMDIndex` exposes bidirectional/FMD search through `SA_Range`. Normal extension delegates to BWA's `bwt_extend()`. The singleton routines are optimized for the case where the current search interval contains exactly one suffix-array row.

The important point is that the singleton optimization does **not** need to reconstruct all four child intervals. Once the requested base is known to occur in the single relevant BWT row, the extension has exactly one row of output.

The companion interval therefore does **not** need the general `child_start()` calculation used for a non-singleton interval.

---

## FMD interval representation

A bidirectional `SA_Range` contains two coupled intervals:

```text
primary interval
companion interval
```

For the BWA FMD representation these are the two sides of BWA's FMD interval. `bwt_extend()` updates both sides together.

The wrapper uses BWA's convention:

```text
is_back = 1  -> logical left extension
is_back = 0  -> logical right extension
```

Thus:

```text
extend_left(c)  -> bwt_extend(..., is_back = 1)
extend_right(c) -> bwt_extend(..., is_back = 0)
```

---

## General extension in `bwt_extend()`

BWA calculates all four possible child intervals. The relevant part is:

```cpp
ok[3].x[is_back] =
    ik->x[is_back] + contains_dollar;

ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
```

Equivalently, the start of child `c` is:

```text
old_start
+ contains_dollar
+ sum(child_size[j] for j > c)
```

This is the correct general rule.

---

## Why the singleton case simplifies

For a singleton, the current interval has:

```text
x[2] = 1
```

There is exactly one current BWT row. For a requested base `c`, that row either contains `c` or it does not.

Therefore:

```text
extension fails -> child_size[c] = 0
extension succeeds -> child_size[c] = 1
```

In the successful case, every other child has size zero:

```text
child_size[j] = 0  for j != c
child_size[c] = 1
```

The cumulative term in the general boundary formula is therefore:

```text
sum(child_size[j] for j > c) = 0
```

because the only non-zero child is the selected child itself.

So the paired boundary reduces to:

```text
new_paired_start = old_paired_start + contains_dollar
```

There is no need for `singleton_extension_sizes()` or `child_start()`.

---

## The `$` row

BWA's FMD index has a special row at:

```text
bwt->primary
```

This row represents the `$` position and is not stored as an ordinary nucleotide in the packed BWT.

The optimized singleton code therefore checks the relevant row explicitly before calling `bwt_B0()`:

```cpp
if (primary == bwt_->primary)
    return false;
```

for left extension, and similarly for the companion side during right extension.

This has an important consequence:

> A successful nucleotide singleton extension cannot have its relevant singleton row equal to `bwt->primary`.

Therefore, on a successful singleton extension:

```text
contains_dollar = 0
```

and the paired interval start is simply unchanged.

The optimized implementation therefore uses the even simpler successful-path form directly: the paired start is unchanged.

---

## Left singleton extension

For:

```cpp
extend_left_singleton(range, c)
```

the BWT operation is performed on the primary interval.

Let:

```text
primary   = p.l
companion = q.l
```

The implementation first checks the single BWT row:

```cpp
if (primary == bwt_->primary)
    return false;

const bwtint_t packed_primary =
    primary - (primary > bwt_->primary);

const uint8_t bwt_c =
    static_cast<uint8_t>(bwt_B0(bwt_, packed_primary));

if (bwt_c != c)
    return false;
```

If the requested base is present, the new primary row is obtained by LF:

```cpp
const bwtint_t new_primary =
    bwt_->L2[c] + 1 +
    bwt_occ(
        bwt_,
        primary == 0 ? (bwtint_t)-1 : primary - 1,
        c);
```

The paired interval is a singleton as well, and the general child-boundary equation has no contribution from other child sizes. Thus:

```cpp
const bwtint_t new_companion = companion;
```

The `$` correction is zero on the successful path, because the `$` row was rejected before the character test.

---

## Right singleton extension

For:

```cpp
extend_right_singleton(range, c)
```

the BWT operation is performed on the companion interval.

The implementation checks the corresponding BWT row:

```cpp
if (companion == bwt_->primary)
    return false;

const bwtint_t packed_companion =
    companion - (companion > bwt_->primary);

const uint8_t bwt_c =
    static_cast<uint8_t>(bwt_B0(bwt_, packed_companion));

if (bwt_c != c)
    return false;
```

The new companion row is then:

```cpp
const bwtint_t new_companion =
    bwt_->L2[c] + 1 +
    bwt_occ(
        bwt_,
        companion == 0 ? (bwtint_t)-1 : companion - 1,
        c);
```

The paired primary boundary follows the singleton simplification:

```cpp
const bwtint_t new_primary = primary;
```

The `$` correction is zero on the successful path.

---

## Why `bwt_occ(l - 1)` is required

BWA's `bwt_occ()` uses an inclusive BWT coordinate.

For a singleton row at `l`, the number of occurrences of base `c` in that row is determined by the transition across the boundary before `l`.

Consequently the LF calculation is:

```cpp
bwt_occ(bwt, l - 1, c)
```

not:

```cpp
bwt_occ(bwt, l, c)
```

For `l == 0`, BWA uses the unsigned representation of `-1`:

```cpp
l == 0 ? (bwtint_t)-1 : l - 1
```

This detail is essential for matching BWA's LF semantics.

---

## Packed BWT and the primary row

BWA's packed BWT does not contain an ordinary nucleotide entry for the special `$` row. Logical BWT row `k` therefore maps to packed row:

```text
k                  if k < bwt->primary
k - 1              if k > bwt->primary
```

The singleton implementation uses:

```cpp
const bwtint_t packed_k =
    k - (k > bwt_->primary);
```

after first rejecting `k == bwt_->primary`.

This is why simply calling `bwt_B0(bwt_, k)` is not sufficient.

---

## Final singleton equations

For a successful left singleton extension:

```text
new_primary   = LF(primary, c)
new_companion = companion
```

For a successful right singleton extension:

```text
new_companion = LF(companion, c)
new_primary   = primary
```

The `$`-aware general form contains an additional `contains_dollar` term, but for a **successful nucleotide singleton extension** that term is necessarily zero because the relevant `$` row is rejected before the extension succeeds.

This is the key simplification compared with the general `bwt_extend()` implementation.

---

## Validation strategy

The optimized routines should be checked against the normal BWA implementation:

```cpp
bwt_extend(bwt, &ik, ok, 1);  // left
bwt_extend(bwt, &ik, ok, 0);  // right
```

For every tested singleton and each base `A/C/G/T`, compare:

```text
optimized singleton result
```

with the corresponding `ok[c]` interval returned by `bwt_extend()`.

The existing diagnostic test was designed around exactly this comparison.

---

## Summary

The general implementation needs the child-size accumulation because several child intervals can be non-empty.

The singleton implementation does not:

1. There is only one current BWT row.
2. A requested base either occurs once or not at all.
3. On success, only the requested child has size one.
4. Therefore all `j > c` child sizes are zero.
5. The cumulative child-size term disappears.
6. The only remaining general correction is the `$` boundary adjustment.
7. Since a `$` singleton cannot be extended by a nucleotide, that adjustment is zero on every successful singleton extension.

Thus the optimized singleton operations reduce to one BWT-row character check plus one LF calculation, with the paired interval boundary carried through unchanged on the successful path.
