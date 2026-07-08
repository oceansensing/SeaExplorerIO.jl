using SeaExplorerIO
using Test
using Dates
using CodecZlib: GzipCompressorStream

# ── synthetic file writers ──────────────────────────────────────────────────

function _write_gz(path, content)
    open(path, "w") do io
        gz = GzipCompressorStream(io)
        write(gz, content)
        close(gz)
    end
end

function _gli_file(dir; name = "sea064.99.gli.sub.1", n = 30)
    t0 = DateTime(2022, 10, 20, 12)
    open(joinpath(dir, name), "w") do io
        write(io, "Timestamp;NavState;Heading;Pitch;Roll;Depth;Lat;Lon;" *
                  "BallastPos;LinPos;Voltage;\n")
        for k in 1:n
            write(io, Dates.format(t0 + Second(10k), "dd/mm/yyyy HH:MM:SS") *
                      ";117;128.06;21.32;8.05;$(5.0 + k);7001.296;201.504;" *
                      "277.4;36.2;28.9;\n")
        end
    end
end

# reference dataset for gated acceptance tests
const M38_DIR = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete"
const M38_NAV = joinpath(M38_DIR, "delayed/nav/logs")
const M38_PLD = joinpath(M38_DIR, "delayed/pld1/logs")

@testset "SeaExplorerIO.jl" begin

    @testset "nmea_to_deg" begin
        @test nmea_to_deg(3908.322) ≈ 39 + 8.322 / 60
        @test nmea_to_deg(-6307.042) ≈ -(63 + 7.042 / 60)
        @test nmea_to_deg(7001.296) ≈ 70 + 1.296 / 60
        @test nmea_to_deg(0.0) == 0.0
        @test isnan(nmea_to_deg(NaN))
    end

    @testset "seaexplorer_files: natural sort, stream isolation" begin
        mktempdir() do d
            for n in (1, 2, 10, 100)
                touch(joinpath(d, "sea064.38.gli.sub.$n.gz"))
            end
            touch(joinpath(d, "sea.0.gli.evt.1"))            # must be ignored
            touch(joinpath(d, "sea064.38.pld1.raw.3.gz"))    # different stream
            fs = seaexplorer_files(d, "gli.sub")
            @test length(fs) == 4
            @test [match(r"\.(\d+)\.gz$", f).captures[1] for f in fs] == ["1", "2", "10", "100"]
            @test seaexplorer_files(d, "pld1.raw") |> length == 1
        end
    end

    @testset "missing_segments" begin
        mktempdir() do d
            for n in (1, 2, 5, 7)
                touch(joinpath(d, "sea064.38.gli.sub.$n.gz"))
            end
            @test missing_segments(d, "gli.sub") == [3, 4, 6]
            @test missing_segments(d, "pld1.raw") == Int[]
        end
    end

    @testset "read_gli: multi-file concat, 1970 filter, NMEA coords" begin
        dir = mktempdir()
        hdr = "Timestamp;NavState;Heading;Pitch;Roll;Depth;Lat;Lon;Voltage;\n"
        # File 2 written first on disk but named .2 — must sort after .1.
        _write_gz(joinpath(dir, "sea064.48.gli.sub.2.gz"),
            hdr *
            "12/11/2023 18:10:00;100;180.00;-20.50;1.00;55.0;7040.100;-800.500;28.0;\n")
        _write_gz(joinpath(dir, "sea064.48.gli.sub.1.gz"),
            hdr *
            "01/01/1970 00:00:10;97;0.00;0.00;0.00;0.0;0.000;0.000;29.0;\n" *   # bench row
            "12/11/2023 18:00:00;105;170.00;25.30;-2.00;5.0;7040.000;-800.400;28.5;\n")

        t = read_gli(dir)
        @test length(t) == 2                              # 1970 row dropped
        @test issorted(t.time)
        @test t.time[1] == DateTime(2023, 11, 12, 18, 0, 0)
        @test t["Pitch"] == [25.30, -20.50]
        @test t["Lat"][1] ≈ 70 + 40.0 / 60                # NMEA → degrees
        @test t["Lon"][2] ≈ -(8 + 0.5 / 60)
        @test length(t.sources) == 2

        # Raw NMEA preserved when conversion is off
        traw = read_gli(dir; convert_coords = false)
        @test traw["Lat"][1] == 7040.000

        # epoch_min = nothing keeps the bench row
        tall = read_gli(dir; epoch_min = nothing)
        @test length(tall) == 3

        # Column subset
        tsub = read_gli(dir; columns = ["Pitch"])
        @test collect(keys(tsub)) == ["Pitch"]
    end

    @testset "read_pld: column selection, sparse rows, ms timestamps" begin
        dir = mktempdir()
        hdr = "PLD_REALTIMECLOCK;NAV_LATITUDE;NAV_LONGITUDE;LEGATO_TEMPERATURE;MR1000G-RDL_EPS1;\n"
        _write_gz(joinpath(dir, "sea064.48.pld1.raw.1.gz"),
            hdr *
            "12/11/2023 18:00:00.100;7040.000;-800.400;3.512;;\n" *
            "12/11/2023 18:00:01.100;7040.001;-800.401;;;\n" *          # all-selected empty
            "12/11/2023 18:00:02.100;7040.002;-800.402;3.514;1.2e-09;\n" *
            "12/11/2023 18:00:03.100;7040.003;-800.403;nan;;\n")        # "nan" cell

        cols = ["LEGATO_TEMPERATURE", "MR1000G-RDL_EPS1"]
        t = read_pld(dir, cols)
        @test length(t) == 3                              # empty row skipped
        @test Dates.millisecond(t.time[1]) == 100
        @test t["LEGATO_TEMPERATURE"][1] ≈ 3.512
        @test t["MR1000G-RDL_EPS1"][2] ≈ 1.2e-9
        @test isnan(t["MR1000G-RDL_EPS1"][1])
        @test isnan(t["LEGATO_TEMPERATURE"][3])

        # skip_empty=false keeps every row; nav coords convert
        tall = read_pld(dir, ["NAV_LATITUDE", "LEGATO_TEMPERATURE"];
                        skip_empty = false)
        @test length(tall) == 4
        @test tall["NAV_LATITUDE"][4] ≈ 70 + 40.003 / 60

        # all columns when none specified
        tfull = read_pld(dir; skip_empty = false)
        @test length(keys(tfull)) == 4

        # Unknown column warns (once, aggregated) and degrades to all-NaN
        tmiss = @test_logs (:warn, r"column absent") match_mode = :any begin
            read_pld(dir, ["NOT_A_COLUMN"]; skip_empty = false)
        end
        @test haskey(tmiss, "NOT_A_COLUMN") && all(isnan, tmiss["NOT_A_COLUMN"])
    end

    @testset "read_stream: generic streams (legato.raw, pld1.sub)" begin
        mktempdir() do d
            hdr = "PLD_REALTIMECLOCK;LEGATO_CONDUCTIVITY;LEGATO_TEMPERATURE;\n"
            _write_gz(joinpath(d, "sea064.48.legato.raw.1.gz"),
                hdr * "12/11/2023 18:00:00.100;38.01;3.512;\n" *
                      "12/11/2023 18:00:01.100;38.02;3.513;\n")
            # a pld1.sub file in the same dir must not be picked up
            _write_gz(joinpath(d, "sea064.48.pld1.sub.1.gz"),
                hdr * "12/11/2023 18:00:00.100;99.0;99.0;\n")
            t = read_stream(d, "legato.raw")
            @test length(t) == 2
            @test t["LEGATO_CONDUCTIVITY"][2] ≈ 38.02
            tsub = read_stream(d, "pld1.sub")
            @test length(tsub) == 1 && tsub["LEGATO_TEMPERATURE"][1] == 99.0
        end
    end

    @testset "mixed headers across files stay row-aligned" begin
        mktempdir() do d
            _write_gz(joinpath(d, "sea064.48.gli.sub.1.gz"),
                "Timestamp;Pitch;Roll;\n" *
                "12/11/2023 18:00:00;10.0;1.0;\n")
            _write_gz(joinpath(d, "sea064.48.gli.sub.2.gz"),
                "Timestamp;Pitch;Voltage;\n" *                   # Roll gone, Voltage new
                "12/11/2023 18:10:00;-11.0;28.5;\n")
            t = read_gli(d)
            @test t["Pitch"] == [10.0, -11.0]
            @test t["Roll"][1] == 1.0 && isnan(t["Roll"][2])
            @test isnan(t["Voltage"][1]) && t["Voltage"][2] == 28.5
        end
    end

    @testset "robustness: missing paths and corrupt files" begin
        # missing directory → actionable error, not a raw open() failure
        @test_throws ErrorException read_gli(joinpath(mktempdir(), "nope"))
        err = try read_gli(joinpath(mktempdir(), "nope")); catch e; e end
        @test occursin("log path not found", sprint(showerror, err))

        # directory with no matching files
        d = mktempdir()
        touch(joinpath(d, "unrelated.txt"))
        err = try read_gli(d); catch e; e end
        @test occursin("no files matching", sprint(showerror, err))

        # a corrupt .gz among good files is skipped; good rows survive;
        # a missing segment number is reported
        d = mktempdir()
        _gli_file(d)
        write(joinpath(d, "sea064.99.gli.sub.3.gz"), UInt8[0x1f, 0x8b, 0xde, 0xad])
        t = @test_logs (:warn, r"missing segment") (:warn, r"unreadable log file skipped") (:warn, r"1 of 2") match_mode = :any begin
            read_gli(d; columns = ["Pitch"])
        end
        @test length(t) == 30 && all(isfinite, t["Pitch"])

        # every file corrupt → hard error
        d = mktempdir()
        write(joinpath(d, "sea064.99.gli.sub.1.gz"), UInt8[0x1f, 0x8b, 0x00])
        err = try read_gli(d); catch e; e end
        @test occursin("unreadable", sprint(showerror, err))
    end

    @testset "robustness: absent columns degrade to NaN" begin
        d = mktempdir()
        _gli_file(d)                                 # has no "Altitude" column
        t = @test_logs (:warn, r"column absent") match_mode = :any begin
            read_gli(d; columns = ["Pitch", "Altitude"])
        end
        @test haskey(t, "Altitude") && all(isnan, t["Altitude"])
        @test all(isfinite, t["Pitch"])
    end

    @testset "Aqua quality assurance" begin
        import Aqua
        Aqua.test_all(SeaExplorerIO)
    end

    if isdir(M38_NAV) && isdir(M38_PLD)
        @testset "M38 acceptance: full-mission nav + payload subset" begin
            nav = read_gli(M38_NAV)
            @test length(nav) > 100_000
            @test issorted(nav.time)
            @test all(v -> isnan(v) || 60 < v < 75, nav["Lat"])   # Lofoten basin, degrees
            @test isempty(missing_segments(M38_NAV, "gli.sub"))
            files = seaexplorer_files(M38_PLD, "pld1.raw")[1:2]
            pld = read_pld(files, ["LEGATO_TEMPERATURE", "AD2CP_HEADING"])
            @test length(pld) > 1000
            temps = filter(isfinite, pld["LEGATO_TEMPERATURE"])
            @test !isempty(temps) && all(-2 .< temps .< 20)
        end
    else
        @info "M38 reference data not found — skipping acceptance tests"
    end

end
