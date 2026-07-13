"""
    SeaExplorerIO

Readers for Alseamar SeaExplorer glider log files: navigation (`*.gli.sub.N[.gz]`)
and payload (`*.pld1.raw.N[.gz]`, `*.pld1.sub.N[.gz]`, `*.legato.raw.N[.gz]`, …)
streams, as logged by the glider and its payload computer — plus the concatenated
`*.<stream>.all.csv` exports served by Alseamar's GLIMPSE command-and-control
server (leading `YO_NUMBER` column, ±9999 fill sentinels, extra derived columns).
Data downloaded several ways can be read together and deduplicated with
[`merge_tables`](@ref) or by passing a vector of directories to the readers.

This package is the file layer only — discovery, gap detection, parsing,
timestamp/coordinate normalization, and robustness to corrupt or missing files.
Sensor physics lives downstream (GliderADCP.jl for ADCP work, GliderTurbulence.jl for
microstructure), where thin wrappers adapt the [`GliderTable`](@ref) returned here.

File conventions (verified on SeaExplorer missions sea064 M37/M38/M48/M59):

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

export GliderTable, read_stream, read_gli, read_pld, merge_tables,
       seaexplorer_files, glimpse_files, missing_segments, nmea_to_deg

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

# Default fill sentinels: SeaExplorer payload logs and GLIMPSE exports write
# ±9999 where an instrument reported nothing (e.g. every AD2CP_* cell while the
# ADCP is off). Disable with `sentinels = nothing`.
const DEFAULT_SENTINELS = (9999.0, -9999.0)

_parse_cell(s::AbstractString, sentinels) = begin
    isempty(s) && return NaN
    v = tryparse(Float64, s)
    v === nothing && return NaN      # "nan", garbage → NaN
    sentinels !== nothing && v in sentinels && return NaN
    return v
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
    glimpse_files(dir, stream) -> Vector{String}

List GLIMPSE server exports of one stream in `dir`, in merge-priority order:
the concatenated whole-mission `<glider>.<mission>.<stream>.all.csv` first,
then per-cycle `<glider>.<mission>.<stream>.<NNN>.csv` exports numerically.
Both flavors carry a leading `YO_NUMBER` column and use ±9999 fill sentinels;
[`read_stream`](@ref) handles them automatically and treats each export as its
own merge source, so overlapping downloads (a cycle also covered by the
`.all.csv`, re-downloaded cycles, decimated exports) deduplicate instead of
doubling rows.
"""
function glimpse_files(dir::AbstractString, stream::AbstractString)
    esc = replace(stream, "." => "\\.")
    pall = Regex("\\." * esc * "\\.all\\.csv\$", "i")
    pcyc = Regex("\\." * esc * "\\.(\\d+)\\.csv\$", "i")
    alls = String[]
    cycs = Tuple{Int,String}[]
    for f in readdir(dir)
        if occursin(pall, f)
            push!(alls, joinpath(dir, f))
        else
            m = match(pcyc, f)
            m === nothing || push!(cycs, (parse(Int, m.captures[1]), joinpath(dir, f)))
        end
    end
    return vcat(sort(alls), last.(sort(cycs)))
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
                         epoch_min::Union{Nothing, DateTime},
                         sentinels)
    io = _open_log(path)
    try
        header = readline(io)
        names = String.(split(strip(header, [';', ' ', '\r']), ';'))
        # The timestamp is column 1 in glider/payload logs; GLIMPSE `.all.csv`
        # exports prepend a YO_NUMBER column, pushing it to column 2 (YO_NUMBER
        # itself is kept as an ordinary data column). Map each selected column
        # NAME to its cell index in THIS file (first occurrence wins on
        # duplicates). Columns are kept row-aligned globally: a column new to
        # this file is back-filled with NaN for all prior rows, and columns
        # absent from this file's header get NaN for its rows — mixed headers
        # across files must never shift values onto wrong timestamps.
        glimpse = length(names) >= 2 && names[1] == "YO_NUMBER"
        tsi = glimpse ? 2 : 1
        src = OrderedDict{String, Int}()
        for (i, n) in enumerate(names)
            i == tsi && continue
            (want === nothing || n in want) || continue
            haskey(src, n) && continue                # duplicate header name
            src[n] = i
            haskey(cols, n) || (cols[n] = fill(NaN, length(times)))
        end
        sel = collect(values(src))
        # "this cell carries data": non-empty, and for GLIMPSE exports (which
        # write ±9999 fills instead of empty cells) also non-sentinel
        carries(parts, i) = i <= length(parts) && !isempty(parts[i]) &&
            (!glimpse || !isnan(_parse_cell(parts[i], sentinels)))

        for line in eachline(io)
            isempty(line) && continue
            parts = split(rstrip(line, '\r'), ';')
            length(parts) <= tsi && continue
            t = _parse_time(parts[tsi])
            t === nothing && continue                 # unparseable row
            epoch_min !== nothing && t < epoch_min && continue   # clock-not-set bench rows

            if skip_empty
                any(i -> carries(parts, i), sel) || continue
            end
            push!(times, t)
            for (name, vec) in cols
                i = get(src, name, 0)
                push!(vec, (i > 0 && i <= length(parts)) ?
                           _parse_cell(parts[i], sentinels) : NaN)
            end
        end
    finally
        close(io)
    end
end

# Read one ordered file list (one source) into a GliderTable.
function _read_source(files::Vector{String}, stream::AbstractString;
                      columns::Union{Nothing, Vector{String}},
                      skip_empty::Bool,
                      epoch_min::Union{Nothing, DateTime},
                      convert_coords::Bool,
                      sentinels,
                      ensure_columns::Bool = true)
    times = DateTime[]
    cols  = OrderedDict{String, Vector{Float64}}()
    nbad = 0
    for f in files
        try
            _read_log_file!(times, cols, f;
                            want = columns, skip_empty, epoch_min, sentinels)
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
    # files existed but yielded nothing (e.g. a zero-byte GLIMPSE export from a
    # failed server download) — say so rather than silently contributing no data
    isempty(times) &&
        @warn "no rows parsed from $(length(files)) file(s) of stream `$stream`" first_file = basename(first(files))

    # Guarantee every requested column exists (all-NaN when the logs never
    # carried it) so downstream indexing degrades instead of KeyError-ing,
    # and say so once rather than once per file. (Deferred to after the merge
    # for multi-source reads — another source may carry the column.)
    if ensure_columns && columns !== nothing
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

"""
    merge_tables(tables::GliderTable...) -> GliderTable

Merge log tables from different download routes (glider computer, GLIMPSE
server) or different resolutions (`pld1.raw`, `pld1.sub`) into one: the union
of all rows and all columns, deduplicated by exact timestamp.

Earlier tables take precedence — list the highest-resolution / most trusted
source first. At a duplicate timestamp the kept row is completed column-wise:
cells that are NaN in the higher-priority row are filled from the
lower-priority one (so e.g. GLIMPSE-computed extra columns attach to
full-resolution rows), and a finite value is never overwritten. Rows whose
timestamp appears in no earlier table are appended whole. Duplicate timestamps
*within* one table are preserved as-is.
"""
function merge_tables(tables::GliderTable...)
    isempty(tables) && error("merge_tables: no tables given")
    length(tables) == 1 && return tables[1]
    allcols = String[]
    for t in tables, k in keys(t)
        k in allcols || push!(allcols, k)
    end
    times = DateTime[]
    cols = OrderedDict{String, Vector{Float64}}(k => Float64[] for k in allcols)
    index = Dict{DateTime, Int}()      # first occurrence, earlier tables only
    ndup = 0
    for t in tables
        added = Dict{DateTime, Int}()
        for i in eachindex(t.time)
            ts = t.time[i]
            j = get(index, ts, 0)
            if j > 0                   # duplicate of an earlier source: coalesce
                ndup += 1
                for k in keys(t)
                    v = t[k][i]
                    isnan(cols[k][j]) && !isnan(v) && (cols[k][j] = v)
                end
            else
                push!(times, ts)
                for k in allcols
                    push!(cols[k], haskey(t, k) ? t[k][i] : NaN)
                end
                haskey(added, ts) || (added[ts] = length(times))
            end
        end
        merge!(index, added)
    end
    if !issorted(times)
        p = sortperm(times)
        times = times[p]
        for (k, v) in cols
            cols[k] = v[p]
        end
    end
    sources = reduce(vcat, (t.sources for t in tables))
    ndup > 0 &&
        @info "merge_tables: $(length(times)) rows from $(length(tables)) sources ($ndup duplicate timestamps coalesced)"
    return GliderTable(times, cols, sources)
end

# The merge sources for one (directory, stream): the sequential segment logs as
# one source, then each GLIMPSE export (.all.csv, then per-cycle .NNN.csv) as its
# OWN source — so any overlap between exports (a cycle also covered by .all.csv,
# re-downloaded or decimated exports) deduplicates at merge time instead of
# doubling rows. Within-source duplicate timestamps stay preserved as always.
function _dir_sources(dir::AbstractString, stream::AbstractString)
    srcs = Vector{String}[]
    segs = seaexplorer_files(dir, stream)
    if !isempty(segs)
        push!(srcs, segs)
        miss = missing_segments(dir, stream)
        isempty(miss) ||
            @warn "SeaExplorer $stream: missing segment numbers" dir = dir missing = miss
    end
    for f in glimpse_files(dir, stream)
        push!(srcs, [f])
    end
    return srcs
end

function _read_logs(src, stream;
                    columns::Union{Nothing, Vector{String}},
                    skip_empty::Bool,
                    epoch_min::Union{Nothing, DateTime},
                    convert_coords::Bool,
                    sentinels)
    streams = stream isa AbstractString ? [String(stream)] : String.(stream)
    kw = (; columns, skip_empty, epoch_min, convert_coords, sentinels)

    src isa AbstractString && !ispath(src) &&
        error("log path not found: $src")

    # explicit single file or file list → one source, first stream label
    if src isa AbstractString && !isdir(src)
        return _read_source([String(src)], streams[1]; kw...)
    end
    if !(src isa AbstractString)
        isempty(src) && error("no files matching stream `$(streams[1])`: empty source list")
        ndirs = count(p -> isdir(String(p)), src)
        if ndirs == 0                          # explicit file list
            return _read_source(String.(src), streams[1]; kw...)
        elseif ndirs < length(src)
            bad = [String(p) for p in src if !isdir(String(p))]
            error("log sources mix directories with files or missing paths: $(join(bad, ", "))")
        end
    end

    # directory source(s) × stream priority list → read each source, merge by
    # priority: streams are ranked (raw before sub), directories break ties in
    # the given order, and within a directory segment logs outrank the GLIMPSE
    # .all.csv, which outranks per-cycle exports. Row union means decimated or
    # fragmentary sources cost nothing: every distinct timestamp from any source
    # survives, so the densest available data wins by construction.
    dirs = src isa AbstractString ? [String(src)] : String.(src)
    tables = GliderTable[]
    nbadsrc = 0
    for st in streams, d in dirs
        for files in _dir_sources(d, st)
            t = try
                _read_source(files, st; kw..., ensure_columns = false)
            catch err
                nbadsrc += 1
                @warn "skipping unreadable source" stream = st nfiles = length(files) error = sprint(showerror, err)
                continue
            end
            push!(tables, t)
        end
    end
    if isempty(tables)
        nbadsrc > 0 &&
            error("all sources of stream `$(join(streams, "`/`"))` in $(join(dirs, ", ")) are unreadable")
        error("no files matching stream `$(join(streams, "`/`"))` in $(join(dirs, ", "))")
    end
    out = merge_tables(tables...)

    if columns !== nothing
        for w in columns
            if !haskey(out, w)
                out.cols[w] = fill(NaN, length(out))
                @warn "column absent from all sources" column = w
            end
        end
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# Public readers
# ─────────────────────────────────────────────────────────────────────────────

const _StreamSpec = Union{AbstractString, AbstractVector{<:AbstractString}}

"""
    read_stream(src, stream; columns = nothing, skip_empty = false,
                epoch_min = DateTime(2000), convert_coords = true,
                sentinels = $(DEFAULT_SENTINELS)) -> GliderTable

Read a SeaExplorer stream (`"gli.sub"`, `"pld1.raw"`, `"pld1.sub"`,
`"legato.raw"`, …) — segment-numbered log files, GLIMPSE-server whole-mission
`.all.csv` exports, and GLIMPSE per-cycle `.NNN.csv` exports are all
recognized. `src` is a directory (with a warning listing any
[`missing_segments`](@ref)), a single file, an explicit vector of file paths,
**or a vector of directories** — any mixture of download routes and export
flavors, possibly decimated or gap-ridden, is combined with
[`merge_tables`](@ref) into one table holding **every distinct timestamp from
any source** (so the densest available data wins by construction), duplicates
deduplicated by priority: streams in the order given (e.g.
`["pld1.raw", "pld1.sub"]` — raw wins), then directories in the order given,
then within a directory segment logs > `.all.csv` > per-cycle exports. At a
shared timestamp the kept row is completed column-wise from lower-priority
sources (GLIMPSE-only derived columns attach to full-resolution rows).

- `columns` selects a subset by header name (nothing = all) — memory stays
  proportional to the selection, which matters for `pld1.raw`
  (10⁶–10⁷ rows × ~60 columns per mission).
- `skip_empty` drops rows where every selected cell is empty (payload
  instruments report asynchronously; in GLIMPSE exports, sentinel-filled
  counts as empty).
- Rows stamped before `epoch_min` are dropped (the glider logs bench rows
  stamped 1970-01-01 before its clock is set); pass `epoch_min = nothing`
  to keep everything.
- Known coordinate columns ($(join(NMEA_COORD_COLS, ", "))) are converted from
  NMEA DDMM.mmm to signed decimal degrees unless `convert_coords = false`.
- Cells equal to a `sentinels` value parse to NaN (SeaExplorer logs and
  GLIMPSE exports use ±9999 as instrument-off fills); `sentinels = nothing`
  disables this.

Corrupt or unreadable files are skipped with a warning; it is an error only
when no file of the stream can be read.
"""
read_stream(src, stream::_StreamSpec;
            columns::Union{Nothing, Vector{String}} = nothing,
            skip_empty::Bool = false,
            epoch_min::Union{Nothing, DateTime} = DateTime(2000),
            convert_coords::Bool = true,
            sentinels = DEFAULT_SENTINELS) =
    _read_logs(src, stream; columns, skip_empty, epoch_min, convert_coords, sentinels)

"""
    read_gli(src; stream = "gli.sub", columns = nothing,
             epoch_min = DateTime(2000), convert_coords = true,
             sentinels = $(DEFAULT_SENTINELS)) -> GliderTable

Read SeaExplorer navigation logs: segment files (`*.gli.sub.N[.gz]`) and/or
GLIMPSE-server exports (`*.gli.sub.all.csv`). `src` is a directory (all
matching files, numerically sorted), a single file, a vector of paths, or a
vector of directories — e.g.
`read_gli(["mission/delayed/nav/logs", "mission/glimpse"])` loads both
download routes and deduplicates by timestamp, earlier directories winning
(see [`merge_tables`](@ref)). `columns` selects a subset by header name
(default: all — the nav table is small). `Lat`/`Lon` are converted from NMEA
DDMM.mmm to decimal degrees unless `convert_coords = false`. Rows stamped
before `epoch_min` (glider clock not yet set) are dropped.
"""
read_gli(src; stream::_StreamSpec = "gli.sub",
         columns::Union{Nothing, Vector{String}} = nothing,
         epoch_min::Union{Nothing, DateTime} = DateTime(2000),
         convert_coords::Bool = true,
         sentinels = DEFAULT_SENTINELS) =
    _read_logs(src, stream; columns, skip_empty = false, epoch_min, convert_coords, sentinels)

"""
    read_pld(src, columns = nothing; stream = ["pld1.raw", "pld1.sub"],
             skip_empty = true, epoch_min = DateTime(2000),
             convert_coords = true, sentinels = $(DEFAULT_SENTINELS)) -> GliderTable

Read SeaExplorer payload logs, keeping only `columns` (header names, e.g.
`["LEGATO_TEMPERATURE", "LEGATO_SALINITY"]`; nothing = all columns — mind the
memory on `pld1.raw`).

`stream` is a priority-ranked list: by default the full-resolution `pld1.raw`
segments are read first and the telemetered `pld1.sub` rows (segment files
and/or GLIMPSE `.all.csv` exports) then fill in only what raw is missing —
all data is loaded, duplicate timestamps keep the highest-resolution values,
and GLIMPSE-only derived columns attach to the merged rows. `src` may be one
directory, a vector of directories (glider-computer + GLIMPSE downloads), a
single file, or an explicit file list.

Because payload instruments report asynchronously, rows where every selected
column is empty (or sentinel-filled, for GLIMPSE exports) are skipped unless
`skip_empty = false`. `NAV_LATITUDE`/`NAV_LONGITUDE` are converted to decimal
degrees unless `convert_coords = false`.
"""
read_pld(src, columns::Union{Nothing, Vector{String}} = nothing;
         stream::_StreamSpec = ["pld1.raw", "pld1.sub"], skip_empty::Bool = true,
         epoch_min::Union{Nothing, DateTime} = DateTime(2000),
         convert_coords::Bool = true,
         sentinels = DEFAULT_SENTINELS) =
    _read_logs(src, stream; columns, skip_empty, epoch_min, convert_coords, sentinels)

end # module
