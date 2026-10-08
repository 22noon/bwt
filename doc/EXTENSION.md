# FM-index Extension in `SA_Range`

## Representation

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

## Extension

The actual extension is performed by BWA's:

```cpp
bwt_extend(bwt_, &ik, ok, is_back);
```

where `ik` is the current FMD interval and `ok[4]` contains the four possible extensions.

The wrapper uses:

```cpp
is_back = 1    → extend P to the left:  cP
is_back = 0    → extend P to the right: Pc
```

For example:

```cpp
SA_Range r2 = index.extend_right(r, c);
```

calls:

```cpp
bwt_extend(bwt_, &ik, ok, 0);
```

and selects `ok[c]`.

Similarly:

```cpp
SA_Range r2 = index.extend_left(r, c);
```

calls:

```cpp
bwt_extend(bwt_, &ik, ok, 1);
```

and selects `ok[c]`.

## Why both intervals are obtained from one operation

For a right extension:

```text
P → Pc
```

the companion pattern changes simultaneously:

```text
RC(P) → RC(Pc)
     = RC(c)RC(P)
```

BWA's FMD representation maintains this paired interval as part of the same extension operation. The wrapper therefore does **not** perform an independent FM search on the companion side.

For a left extension:

```text
P → cP
```

the companion becomes:

```text
RC(P) → RC(cP)
     = RC(P)RC(c)
```

Again, `bwt_extend()` maintains both coordinates.

## Branching

The `extend_left_all()` and `extend_right_all()` operations call `bwt_extend()` once and retain all four results:

```cpp
bwtintv_t ok[4]{};
bwt_extend(bwt_, &ik, ok, is_back);
```

This is the appropriate primitive when the search branches because all four possible nucleotide extensions are produced together.

## Singleton path

For a singleton interval (`size == 1`), constructing all four `ok[]` intervals is unnecessary.

The planned singleton fast path will instead:

1. inspect the BWT character at the unique primary row;
2. reject immediately if it is not the requested base;
3. perform the corresponding LF calculation;
4. derive the paired FMD coordinate without independently extending the companion interval.

This provides a cheaper primitive for the common case where the search follows a unique path.

