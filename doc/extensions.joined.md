# FM-index Extension in `SA_Range` and Singleton Extension in the BWA FMD Index

## 1. FM-index Extension in `SA_Range`

### Representation

`SA_Range` represents an FMD/bidirectional search state using two intervals:

```text
primary     = SA interval for the current pattern P
companion   = SA interval for RC(P)
```

For an interval `[l,r)`, the size is `r-l`.

Internally, the BWA backend represents the same state as:

```cpp
bwtintv_t {
    x[0] = primary.l,
    x[1] = companion.l,
    x[2] = interval size
}
```

Thus conversion between `SA_Range` and BWA's `bwtintv_t` does not require searching.

---

## Initialisation

For a single nucleotide `c`:

```cpp
bwtintv_t v{};
bwt_set_intv(bwt_, c, v);
```

BWA constructs the initial FMD interval. The wrapper converts:

```text
v.x[0]       → primary interval
v.x[1]       → companion interval
v.x[2]       → size
```

into an `SA_Range`.

---

## Extension

The actual extension is performed by BWA's:

```cpp
bwt_extend(bwt_, &ik, ok, is_back);
```

where `ik` is the current FMD interval and `ok[4]` contains the four possible BWA/FMD child intervals.

The public `BwaFMDIndex` API uses biological/logical extension semantics:

```text
extend_left(c)  -> cP
extend_right(c) -> Pc
```

The mapping to BWA's FMD representation is:

```text
logical operation       is_back       BWA child

extend_left(c)           1             ok[c]

extend_right(c)          0             ok[complement(c)]
```

where BWA's nucleotide encoding is:

```text
A=0, C=1, G=2, T=3
```

and therefore:

```cpp
complement(c) = 3 - c;
```

The complement is required because BWA's FMD companion side represents the reverse complement of the primary pattern.

For example, for:

```text
P = AC
```

a logical right extension by `G` produces:

```text
AC -> ACG
```

while the companion pattern changes from:

```text
RC(AC) = GT
```

to:

```text
RC(ACG) = CGT
```

Thus BWA's `is_back=0` operation uses the complementary nucleotide when selecting the corresponding FMD child.

The wrapper hides this representation detail:

```cpp
SA_Range r2 = index.extend_right(r, c);
```

always means "append biological base `c` to the current pattern".

The caller does not need to apply the complement itself.

Similarly:

```cpp
SA_Range r2 = index.extend_left(r, c);
```

means "prepend biological base `c`".

---

## Why both intervals are obtained from one operation

For a right extension:

```text
P -> Pc
```

the companion pattern changes simultaneously:

```text
RC(P) -> RC(Pc)
     = RC(c)RC(P)
```

BWA's FMD representation maintains this paired interval as part of the same extension operation. The wrapper therefore does **not** perform an independent FM search on the companion side.

For a left extension:

```text
P -> cP
```

the companion becomes:

```text
RC(P) -> RC(cP)
     = RC(P)RC(c)
```

Again, `bwt_extend()` maintains both coordinates.

---

## Branching

The `extend_left_all()` and `extend_right_all()` operations call `bwt_extend()` once and retain all four results.

### `extend_left_all()`

```text
out[A] = BWA ok[A]
out[C] = BWA ok[C]
out[G] = BWA ok[G]
out[T] = BWA ok[T]
```

### `extend_right_all()`

```text
out[A] = BWA ok[T]
out[C] = BWA ok[G]
out[G] = BWA ok[C]
out[T] = BWA ok[A]
```

because a logical right-extension base is complemented before selecting the BWA FMD child.

This is the appropriate primitive when the search branches because all four possible nucleotide extensions are produced together.

---

# 2. Singleton Extension in the BWA FMD Index

## Overview

`BwaFMDIndex` exposes bidirectional/FMD search through `SA_Range`. Normal extension delegates to BWA's `bwt_extend()`. The singleton routines are optimized for the case where the current search interval contains exactly one suffix-array row.

The important point is that the singleton optimization does **not** need to reconstruct all four child intervals. Once the requested base is known to occur in the single relevant BWT row, the extension has exactly one row of output.

The companion interval therefore does not need the general child-start calculation used for a non-singleton interval.

---

## FMD interval representation

A bidirectional `SA_Range` contains two coupled intervals:

```text
primary interval
companion interval
```

For the BWA FMD representation:

```text
primary    = interval for P
companion  = interval for RC(P)
```

`bwt_extend()` updates both sides together.

The wrapper uses BWA's convention:

```text
is_back = 1  -> logical left extension
is_back = 0  -> logical right extension
```

However, the nucleotide selected from BWA's `ok[]` array differs between the two directions:

```text
logical left extension by c
    -> ok[c]

logical right extension by c
    -> ok[complement(c)]
```

This distinction is essential because the companion side is reverse-complemented.

---

## General extension in `bwt_extend()`

BWA calculates all four possible child intervals. The relevant boundary calculation is:

```cpp
ok[3].x[is_back] =
    ik->x[is_back] + contains_dollar;

ok[2].x[is_back] =
    ok[3].x[is_back] + ok[3].x[2];

ok[1].x[is_back] =
    ok[2].x[is_back] + ok[2].x[2];

ok[0].x[is_back] =
    ok[1].x[is_back] + ok[1].x[2];
```

Equivalently, the start of child `c` is:

```text
old_start
+ contains_dollar
+ sum(child_size[j] for j > c)
```

For a logical right extension, remember that the BWA child is indexed by `complement(c)`.

---

## Why the singleton case simplifies

For a singleton, the current interval has:

```text
x[2] = 1
```

There is exactly one current BWT row. For a requested logical base `c`, that row either contains the required BWT character or it does not.

Therefore:

```text
extension fails    -> selected child_size = 0
extension succeeds -> selected child_size = 1
```

In the successful case, every other child has size zero:

```text
child_size[j] = 0  for j != selected_child
child_size[selected_child] = 1
```

The cumulative term in the general boundary formula is therefore:

```text
sum(child_size[j] for j > selected_child) = 0
```

because the only non-zero child is the selected child itself.

So the paired boundary reduces to:

```text
new_paired_start = old_paired_start + contains_dollar
```

There is no need for a general child-size calculation.

---

## The `$` row

BWA's FMD index has a special row at:

```cpp
bwt->primary
```

This row represents the `$` position and is not stored as an ordinary nucleotide in the packed BWT.

The optimized singleton code therefore checks the relevant row explicitly before calling `bwt_B0()`.

For left extension:

```cpp
if (primary == bwt_->primary)
    return false;
```

For right extension:

```cpp
if (companion == bwt_->primary)
    return false;
```

A successful nucleotide singleton extension therefore cannot have its relevant singleton row equal to `bwt->primary`.

Consequently:

```text
contains_dollar = 0
```

on every successful singleton nucleotide extension, and the paired interval start is unchanged.

---

# Left singleton extension

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
    static_cast<uint8_t>(
        bwt_B0(bwt_, packed_primary));

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

The successful logical left singleton extension is therefore:

```text
new_primary   = LF(primary, c)
new_companion = companion
```

---

# Right singleton extension

For:

```cpp
extend_right_singleton(range, c)
```

the BWT operation is performed on the companion interval.

The important distinction is that `c` is the **logical biological base**, while the BWT character on the companion side is:

```cpp
const uint8_t fmd_c =
    static_cast<uint8_t>(3 - c);
```

because the companion represents `RC(P)`.

The implementation therefore checks:

```cpp
if (companion == bwt_->primary)
    return false;

const bwtint_t packed_companion =
    companion - (companion > bwt_->primary);

const uint8_t bwt_c =
    static_cast<uint8_t>(
        bwt_B0(bwt_, packed_companion));

if (bwt_c != fmd_c)
    return false;
```

The new companion row is obtained by LF using the FMD character:

```cpp
const bwtint_t new_companion =
    bwt_->L2[fmd_c] + 1 +
    bwt_occ(
        bwt_,
        companion == 0 ? (bwtint_t)-1 : companion - 1,
        fmd_c);
```

The paired primary boundary follows the singleton simplification:

```cpp
const bwtint_t new_primary = primary;
```

Thus the successful logical right singleton extension is:

```text
fmd_c         = complement(c)
new_companion = LF(companion, fmd_c)
new_primary   = primary
```

The caller still supplies only the logical base `c`.

---

## Why `bwt_occ(l - 1)` is required

BWA's `bwt_occ()` uses an inclusive BWT coordinate.

For a singleton row at `l`, the LF calculation must count occurrences strictly before that row. Consequently the calculation is:

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

BWA's packed BWT does not contain an ordinary nucleotide entry for the special `$` row.

Logical BWT row `k` therefore maps to packed row:

```text
k       if k < bwt->primary
k - 1   if k > bwt->primary
```

The singleton implementation uses:

```cpp
const bwtint_t packed_k =
    k - (k > bwt_->primary);
```

after first rejecting:

```cpp
k == bwt_->primary
```

This is why simply calling:

```cpp
bwt_B0(bwt_, k)
```

is not sufficient.

---

# Final singleton equations

For a successful logical left singleton extension:

```text
new_primary   = LF(primary, c)
new_companion = companion
```

For a successful logical right singleton extension:

```text
fmd_c         = complement(c)
new_companion = LF(companion, fmd_c)
new_primary   = primary
```

The `$`-aware general form contains an additional `contains_dollar` term, but for a successful nucleotide singleton extension that term is necessarily zero because the relevant `$` row is rejected before the extension succeeds.

---

# Validation strategy

The optimized routines should be checked against the normal BWA implementation:

```cpp
bwt_extend(bwt, &ik, ok, 1);  // left
bwt_extend(bwt, &ik, ok, 0);  // right
```

For every tested singleton and each logical base `A/C/G/T`, compare:

```text
optimized singleton result
```

with the corresponding BWA child:

```text
logical left extension by c
    -> ok[c]

logical right extension by c
    -> ok[complement(c)]
```

The same mapping must be used by the non-singleton diagnostic tests.

This validates both:

1. the optimized singleton implementation, and
2. the wrapper's logical nucleotide-extension contract.

---

# Summary

The general implementation needs the child-size accumulation because several child intervals can be non-empty.

The singleton implementation does not:

1. There is only one current BWT row.
2. A requested base either occurs once or not at all.
3. On success, only the selected child has size one.
4. Therefore all child sizes greater than the selected child are zero.
5. The cumulative child-size term disappears.
6. The only remaining general correction is the `$` boundary adjustment.
7. Since a `$` singleton cannot be extended by a nucleotide, that adjustment is zero on every successful singleton extension.

The resulting optimized singleton operations reduce to one BWT-row character check plus one LF calculation, with the paired interval boundary carried through unchanged on the successful path.

For BWA's FMD representation, the critical direction-dependent mapping is:

```text
logical operation       BWA child

left extension by c     ok[c]

right extension by c    ok[complement(c)]
```

The `BwaFMDIndex` wrapper hides this representation detail from callers.

