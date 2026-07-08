# SeaExplorerIO.jl

Pure-Julia readers for [Alseamar SeaExplorer](https://www.alseamar-alcen.com/products/underwater-glider/seaexplorer)
glider log files: navigation (`*.gli.sub.N[.gz]`) and payload
(`*.pld1.raw.N[.gz]`, `*.pld1.sub.N[.gz]`, `*.legato.raw.N[.gz]`, …) streams.

This package is deliberately the **file layer only** — discovery, transfer-gap
detection, parsing, timestamp/coordinate normalization, and robustness to
corrupt or missing files. Sensor physics lives downstream:

- [GliderADCP.jl](https://github.com/truedichotomy/GliderADCP.jl) — Nortek
  AD2CP processing to absolute ocean velocities
- [ATOMIXjulia.jl](https://github.com/truedichotomy/ATOMIXjulia.jl) — MicroRider
  microstructure/turbulence processing

Both wrap this package, so a loader bugfix or a new sensor lands once, here.

## Usage

```julia
using SeaExplorerIO

# navigation: all segments in a directory, naturally sorted, time-ordered
nav = read_gli("mission/nav/logs")
nav["Pitch"]                       # Float64 vector, NaN = missing
nav.time                           # DateTime per row
nav["Lat"]                         # NMEA DDMM.mmm auto-converted to degrees

# payload: select columns — memory stays proportional to the selection
pld = read_pld("mission/pld1/logs", ["LEGATO_TEMPERATURE", "LEGATO_SALINITY"])

# any other segment-numbered stream
leg = read_stream("mission/pld1/logs", "legato.raw")

# file bookkeeping
seaexplorer_files("mission/nav/logs", "gli.sub")   # naturally sorted segment files
missing_segments("mission/nav/logs", "gli.sub")    # gaps in the transfer sequence
```

All readers return a `GliderTable` (row timestamps + name → `Vector{Float64}`
columns; `NaN` marks empty/non-numeric cells). Rows stamped before the glider's
clock is set (epoch-1970 bench records) are dropped by default.

### Robustness contract

Real missions are messy; the readers degrade gracefully and loudly:

- corrupt or unreadable segment files are **skipped with a warning** (an error
  is raised only if *no* file of the stream parses)
- missing segment numbers in a transfer sequence are **reported**
- a column absent from every file **warns once** and returns all-NaN instead of
  a `KeyError`; mixed headers across segment files stay row-aligned
- unparseable rows and cells become skipped rows / NaN, never exceptions

## Provenance

Merged from the independently developed loaders of ATOMIXjulia.jl (lean
column-selective parser, schema guarantees) and GliderADCP.jl (generic stream
discovery, transfer-gap detection), keeping the best of both. Validated
against the sea064 M38 NorSE mission (Lofoten Basin, 2022–2023).
