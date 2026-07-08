"""
    SeaExplorerIO

Readers for Alseamar SeaExplorer glider log files: navigation (`*.gli.sub.N[.gz]`)
and payload (`*.pld1.raw.N[.gz]`, `*.pld1.sub.N[.gz]`, `*.legato.raw.N[.gz]`, …)
streams, as logged by the glider and its payload computer.

This package is the file layer only — discovery, gap detection, parsing,
timestamp/coordinate normalization, and robustness to corrupt or missing files.
Sensor physics lives downstream (GliderADCP.jl for ADCP work, ATOMIXjulia.jl for
microstructure), where thin wrappers adapt the [`GliderTable`](@ref) returned here.

File conventions (verified on SeaExplorer missions sea064 M38/M48):

  - semicolon-separated text with a header line and a trailing `;` per row,
    optionally gzipped
  - segment-numbered names `<glider>.<mission>.<stream>.<N>[.gz]` — one file per
    glider yo/transfer; `N` needs natural (numeric) sort
  - nav `Timestamp` = `dd/mm/yyyy HH:MM:SS`; payload `PLD_REALTIMECLOCK` adds `.sss`
  - Lat/Lon columns are NMEA-style DDMM.mmm ([`nmea_to_deg`](@ref) converts)
  - rows stamped before the glider's clock is set carry epoch-1970 dates
    (dropped by default via `epoch_min`)
  - payload instruments report asynchronously, so most payload cells are empty
    on any given row
"""
module SeaExplorerIO

using Dates
using CodecZlib: GzipDecompressorStream
using OrderedCollections: OrderedDict

export GliderTable, read_stream, read_gli, read_pld,
       seaexplorer_files, missing_segments, nmea_to_deg

# ─────────────────────────────────────────────────────────────────────────────
# GliderTable — the neutral column table both downstream packages adapt
# ─────────────────────────────────────────────────────────────────────────────

"""
    GliderTable

Column table of glider log data. `time` is the row timestamp; `cols` maps
column name → Float64 vector (NaN = missing/empty). Index with
`tbl["Pitch"]`; column names via `keys(tbl)`; `length(tbl)` = row count.
`sources` records the files that were read, in concatenation order.

Values are Float64 throughout (including flags like NavState) — this is a
flight/science table, not a bit-fidelity path; non-numeric cells parse to NaN.
"""
struct GliderTable
    time::Vector{DateTime}
    cols::OrderedDict{String, Vector{Float64}}
    sources::Vector{String}
end

Base.getindex(t::GliderTable, name::AbstractString) = t.cols[name]
Base.haskey(t::GliderTable, name::AbstractString)   = haskey(t.cols, name)
Base.keys(t::GliderTable)                           = keys(t.cols)
Base.length(t::GliderTable)                         = length(t.time)

function Base.show(io::IO, t::GliderTable)
    span = isempty(t.time) ? "empty" :
        string(first(t.time), " → ", last(t.time))
    print(io, "GliderTable(", length(t.time), " rows, ",
          length(t.cols), " cols, ", span, ")")
end

# ─────────────────────────────────────────────────────────────────────────────
# Cell/coordinate/timestamp parsing
# ─────────────────────────────────────────────────────────────────────────────

"""
    nmea_to_deg(x) -> Float64

Convert an NMEA-style coordinate DDMM.mmm (degrees × 100 + decimal minutes,
sign carrying the hemisphere) to signed decimal degrees.
`-6307.042 → -(63 + 7.042/60)`. NaN passes through.
"""
function nmea_to_deg(x::Real)
    isnan(x) && return Float64(NaN)
    s = sign(x)
    a = abs(x)
    d = floor(a / 100)
    return s * (d + (a - 100d) / 60)
end

# Column names (lowercase-matched) holding NMEA DDMM.mmm coordinates.
const NMEA_COORD_COLS = ("lat", "lon", "nav_latitude", "nav_longitude")

# Timestamp formats seen across streams (payload adds milliseconds).
const _FMT_SEC = dateformat"dd/mm/yyyy HH:MM:SS"
const _FMT_MS  = dateformat"dd/mm/yyyy HH:MM:SS.sss"

_parse_time(s::AbstractString) =
    tryparse(DateTime, s, length(s) > 19 ? _FMT_MS : _FMT_SEC)

_parse_cell(s::AbstractString) = begin
    isempty(s) && return NaN
    v = tryparse(Float64, s)
    v === nothing ? NaN : v          # "nan", garbage → NaN
end

_open_log(path::AbstractString) = endswith(path, ".gz") ?
    GzipDecompressorStream(open(path, "r")) : open(path, "r")

# ─────────────────────────────────────────────────────────────────────────────
# File discovery and transfer-gap detection
# ─────────────────────────────────────────────────────────────────────────────

_stream_pattern(stream::AbstractString) =
    Regex("\\." * replace(stream, "." => "\\.") * "\\.(\\d+)(\\.gz)?\$")

"""
    seaexplorer_files(dir, stream) -> Vector{String}

List segment files of one SeaExplorer stream (e.g. `"gli.sub"`, `"pld1.raw"`,
`"legato.raw"`, `"ad2cp.raw"`) in `dir`, naturally sorted by segment number.
"""
function seaexplorer_files(dir::AbstractString, stream::AbstractString)
    pat = _stream_pattern(stream)
    hits = Tuple{Int,String}[]
    for f in readdir(dir)
        m = match(pat, f)
        m === nothing || push!(hits, (parse(Int, m.captures[1]), joinpath(dir, f)))
    end
    return last.(sort(hits))
end

"""
    missing_segments(dir, stream) -> Vector{Int}

Segment numbers absent from an otherwise consecutive SeaExplorer stream sequence
(e.g. `sea064.38.gli.sub.<N>` with N = 1…max): the gaps in mission file transfer.
"""
function missing_segments(dir::AbstractString, stream::AbstractString)
    pat = _stream_pattern(stream)
    nums = Int[]
    for f in readdir(dir)
        m = match(pat, f)
        m === nothing || push!(nums, parse(Int, m.captures[1]))
    end
    isempty(nums) && return Int[]
    return setdiff(minimum(nums):maximum(nums), nums)
end

# ─────────────────────────────────────────────────────────────────────────────
# Core reader
# ─────────────────────────────────────────────────────────────────────────────

"""
Parse one semicolon-separated log file into (times, per-column values).
`want` restricts to a column subset (nothing = all columns). Rows where every
selected cell is empty are skipped when `skip_empty` is set.
"""
function _read_log_file!(times::Vector{DateTime},
                         cols::OrderedDict{String, Vector{Float64}},
                         path::AbstractString;
                         want::Union{Nothing, Vector{String}},
                         skip_empty::Bool,
                         epoch_min::Union{Nothing, DateTime})
    io = _open_log(path)
    try
        header = readline(io)
        names = String.(split(strip(header, [';', ' ', '\r']), ';'))
        # Column 1 is the timestamp. Map each selected column NAME to its
        # cell index in THIS file (first occurrence wins on duplicates).
        # Columns are kept row-aligned globally: a column new to this file
        # is back-filled with NaN for all prior rows, and columns absent
        # from this file's header get NaN for its rows — mixed headers
        # across files must never shift values onto wrong timestamps.
        src = OrderedDict{String, Int}()
        for (i, n) in enumerate(names)
            i == 1 && continue
            (want === nothing || n in want) || continue
            haskey(src, n) && continue                # duplicate header name
            src[n] = i
            haskey(cols, n) || (cols[n] = fill(NaN, length(times)))
        end
        sel = collect(values(src))

        for line in eachline(io)
            isempty(line) && continue
            parts = split(rstrip(line, '\r'), ';')
            length(parts) < 2 && continue
            t = _parse_time(parts[1])
            t === nothing && continue                 # unparseable row
            epoch_min !== nothing && t < epoch_min && continue   # clock-not-set bench rows

            if skip_empty
                any(i -> i <= length(parts) && !isempty(parts[i]), sel) || continue
            end
            push!(times, t)
            for (name, vec) in cols
                i = get(src, name, 0)
                push!(vec, (i > 0 && i <= length(parts)) ?
                           _parse_cell(parts[i]) : NaN)
            end
        end
    finally
        close(io)
    end
end

function _read_logs(src, stream::AbstractString;
                    columns::Union{Nothing, Vector{String}},
                    skip_empty::Bool,
                    epoch_min::Union{Nothing, DateTime},
                    convert_coords::Bool)
    src isa AbstractString && !ispath(src) &&
        error("log path not found: $src")
    local files
    if src isa AbstractString && isdir(src)
        files = seaexplorer_files(src, stream)
        isempty(files) &&
            error("no files matching stream `$stream` in $src")
        miss = missing_segments(src, stream)
        isempty(miss) ||
            @warn "SeaExplorer $stream: missing segment numbers" missing = miss
    else
        files = src isa AbstractString ? [String(src)] : String.(src)
        isempty(files) && error("no files matching stream `$stream`: empty file list")
    end

    times = DateTime[]
    cols  = OrderedDict{String, Vector{Float64}}()
    nbad = 0
    for f in files
        try
            _read_log_file!(times, cols, f;
                            want = columns, skip_empty, epoch_min)
        catch err
            # a corrupt/truncated file (e.g. bad gzip) must not take the
            # rest of the mission down with it
            nbad += 1
            @warn "unreadable log file skipped" file=basename(f) error=sprint(showerror, err)
            # drop any rows the partial read may have appended unevenly
            nt = length(times)
            for (k, v) in cols
                length(v) > nt ? resize!(v, nt) :
                length(v) < nt && append!(v, fill(NaN, nt - length(v)))
            end
        end
    end
    nbad == length(files) &&
        error("all $(length(files)) log files of stream `$stream` are unreadable")
    nbad > 0 && @warn "$nbad of $(length(files)) log files skipped as unreadable"

    # Guarantee every requested column exists (all-NaN when the logs never
    # carried it) so downstream indexing degrades instead of KeyError-ing,
    # and say so once rather than once per file.
    if columns !== nothing
        for w in columns
            if !haskey(cols, w)
                cols[w] = fill(NaN, length(times))
                @warn "column absent from all $(length(files)) log files" column=w
            end
        end
    end

    # Files are per-yo and internally ordered, but guarantee global order.
    if !issorted(times)
        p = sortperm(times)
        times = times[p]
        for (k, v) in cols
            cols[k] = v[p]
        end
    end

    if convert_coords
        for (k, v) in cols
            lowercase(k) in NMEA_COORD_COLS && map!(nmea_to_deg, v, v)
        end
    end

    return GliderTable(times, cols, files)
end

# ─────────────────────────────────────────────────────────────────────────────
# Public readers
# ─────────────────────────────────────────────────────────────────────────────

"""
    read_stream(src, stream; columns = nothing, skip_empty = false,
                epoch_min = DateTime(2000), convert_coords = true) -> GliderTable

Read any segment-numbered SeaExplorer stream (`"gli.sub"`, `"pld1.raw"`,
`"pld1.sub"`, `"legato.raw"`, …). `src` is a directory (all matching segment
files, naturally sorted — with a warning listing any [`missing_segments`](@ref)),
a single file, or an explicit vector of paths.

- `columns` selects a subset by header name (nothing = all) — memory stays
  proportional to the selection, which matters for `pld1.raw`
  (10⁶–10⁷ rows × ~60 columns per mission).
- `skip_empty` drops rows where every selected cell is empty (payload
  instruments report asynchronously).
- Rows stamped before `epoch_min` are dropped (the glider logs bench rows
  stamped 1970-01-01 before its clock is set); pass `epoch_min = nothing`
  to keep everything.
- Known coordinate columns ($(join(NMEA_COORD_COLS, ", "))) are converted from
  NMEA DDMM.mmm to signed decimal degrees unless `convert_coords = false`.

Corrupt or unreadable files are skipped with a warning; it is an error only
when no file of the stream can be read.
"""
read_stream(src, stream::AbstractString;
            columns::Union{Nothing, Vector{String}} = nothing,
            skip_empty::Bool = false,
            epoch_min::Union{Nothing, DateTime} = DateTime(2000),
            convert_coords::Bool = true) =
    _read_logs(src, stream; columns, skip_empty, epoch_min, convert_coords)

"""
    read_gli(src; stream = "gli.sub", columns = nothing,
             epoch_min = DateTime(2000), convert_coords = true) -> GliderTable

Read SeaExplorer navigation logs (`*.gli.sub.N[.gz]`). `src` is a directory
(all matching files, numerically sorted), a single file, or a vector of
paths. `columns` selects a subset by header name (default: all — the nav
table is small). `Lat`/`Lon` are converted from NMEA DDMM.mmm to decimal
degrees unless `convert_coords = false`. Rows stamped before `epoch_min`
(glider clock not yet set) are dropped.
"""
read_gli(src; stream::AbstractString = "gli.sub",
         columns::Union{Nothing, Vector{String}} = nothing,
         epoch_min::Union{Nothing, DateTime} = DateTime(2000),
         convert_coords::Bool = true) =
    _read_logs(src, stream; columns, skip_empty = false, epoch_min, convert_coords)

"""
    read_pld(src, columns = nothing; stream = "pld1.raw", skip_empty = true,
             epoch_min = DateTime(2000), convert_coords = true) -> GliderTable

Read SeaExplorer payload logs (`*.pld1.raw.N[.gz]` by default; set `stream`
for `"pld1.sub"`, `"legato.raw"`, …), keeping only `columns` (header names,
e.g. `["LEGATO_TEMPERATURE", "LEGATO_SALINITY"]`; nothing = all columns —
mind the memory on `pld1.raw`). Because payload instruments report
asynchronously, rows where every selected column is empty are skipped unless
`skip_empty = false`. `NAV_LATITUDE`/`NAV_LONGITUDE` are converted to decimal
degrees unless `convert_coords = false`.
"""
read_pld(src, columns::Union{Nothing, Vector{String}} = nothing;
         stream::AbstractString = "pld1.raw", skip_empty::Bool = true,
         epoch_min::Union{Nothing, DateTime} = DateTime(2000),
         convert_coords::Bool = true) =
    _read_logs(src, stream; columns, skip_empty, epoch_min, convert_coords)

end # module
