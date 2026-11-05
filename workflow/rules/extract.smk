import shlex


def _ingest_preextracted(src_dir: Path, dst_dir: Path) -> int:
    """Symlink {peak}.fa/.tsv pairs (e.g. from PARE) into dst_dir and write peak_ids.txt.
    Returns the number of pairs linked."""
    src_dir = Path(src_dir).resolve()
    dst_dir.mkdir(parents=True, exist_ok=True)
    linked: list[str] = []
    for fa in sorted(src_dir.glob("*.fa")):
        peak = fa.stem
        tsv = src_dir / f"{peak}.tsv"
        if not tsv.exists():
            continue
        for src, dst in ((fa, dst_dir / fa.name), (tsv, dst_dir / tsv.name)):
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src)
        linked.append(peak)
    (dst_dir / "peak_ids.txt").write_text("\n".join(linked) + "\n")
    return len(linked)


def _haplotypes_parquet(fastas_dir: Path, out: Path) -> None:
    """Consolidate per-peak fastas (.fa seq + .tsv actuation/coords) into one
    haplotypes.parquet -- the preextracted-mode equivalent of what
    extract_fastas.py writes inline on the union-BED path."""
    import polars as pl

    def read_fa(p: Path) -> dict[str, str]:
        d: dict[str, list[str]] = {}
        name = None
        for ln in p.read_text().splitlines():
            if ln.startswith(">"):
                name = ln[1:].split()[0]
                d[name] = []
            elif name is not None:
                d[name].append(ln.strip())
        return {k: "".join(v) for k, v in d.items()}

    cols = [
        "peak_id",
        "sample",
        "seq",
        "actuation",
        "coverage",
        "fire_coverage",
        "asm_chr",
        "asm_start",
        "asm_end",
    ]
    rows: dict[str, list] = {c: [] for c in cols}
    for fa in sorted(fastas_dir.glob("*.fa")):
        peak = fa.stem
        tsv = fastas_dir / f"{peak}.tsv"
        if not tsv.exists():
            continue
        seqs = read_fa(fa)
        lines = tsv.read_text().splitlines()
        ix = {h: i for i, h in enumerate(lines[0].split("\t"))}
        for ln in lines[1:]:
            f = ln.split("\t")
            s = f[ix["sample"]]
            for k, v in zip(
                cols,
                (
                    peak,
                    s,
                    seqs.get(s, ""),
                    f[ix["actuation"]],
                    f[ix["coverage"]],
                    f[ix["fire_coverage"]],
                    f[ix["asm_chr"]],
                    f[ix["asm_start"]],
                    f[ix["asm_end"]],
                ),
            ):
                rows[k].append(v)
    out.parent.mkdir(parents=True, exist_ok=True)
    (
        pl.DataFrame(rows, schema={c: pl.Utf8 for c in cols})
        .with_columns(
            pl.col("actuation").cast(pl.Float64, strict=False),
            pl.col("coverage").cast(pl.Float64, strict=False),
            pl.col("fire_coverage").cast(pl.Float64, strict=False),
        )
        .write_parquet(out)
    )


rule extract_fastas:
    """Stream the union BED into per-peak FASTA/TSV, haplotypes.parquet and the survivor
manifest, or symlink a preextracted_dir (e.g. PARE paralog groups). Not a checkpoint."""
    input:
        src=_extract_fastas_input,
        peak_ids=_peak_ids_input,
        highlight=_highlight_input,
    output:
        # Every requested peak is declared so the DAG is static. A peak that yields
        # nothing leaves a declared output missing; extract_fastas.py exits listing
        # those ids. fastas/ is temp; haplotypes.parquet is the persisted copy.
        fa=temp(expand(RESULTS + "/fastas/{peak_id}.fa", peak_id=_run_peak_ids())),
        tsv=temp(expand(RESULTS + "/fastas/{peak_id}.tsv", peak_id=_run_peak_ids())),
        manifest=PEAKS,
        haplotypes=HAP,
    threads: config.get("bgzip_threads", 16)
    resources:
        mem_mb=32768,
        runtime=480,
    params:
        script=script_path("extract_fastas.py"),
        # run: blocks don't track config reads; params make config changes re-run this rule.
        preextracted=config.get("preextracted_dir"),
        min_seqs=config.get("min_seqs", 2),
        primary_only=config.get("primary_only", True),
        haplotype_distinct=config.get("haplotype_distinct"),
    run:
        if params.preextracted:
            n = _ingest_preextracted(Path(input.src), Path(RESULTS + "/fastas"))
            _haplotypes_parquet(Path(RESULTS + "/fastas"), Path(output.haplotypes))
            print(f"linked {n} pre-extracted FASTA/TSV pairs from {input.src}")
        else:
            # The union BED is asm-coord-sorted; group rows by consensus_peak_id
            # (col 5) for the streaming extractor. `sed -u 1q` keeps the header first.
            parts = [
                f"bgzip -cd -@ {threads} {shlex.quote(str(input.src))}",
                "| (sed -u 1q; sort -t$'\\t' -k5,5 -S 50%)",
                f"| python {shlex.quote(str(params.script))}",
                f"--outdir {RESULTS}/fastas",
                f"--min-seqs {params.min_seqs}",
                f"--haplotypes-parquet {output.haplotypes}",
            ]
            if params.primary_only:
                parts.append("--primary-only")
            if input.peak_ids:
                parts.append(f"--peak-ids {shlex.quote(str(input.peak_ids[0]))}")
            if params.haplotype_distinct:
                parts.append(
                    f"--haplotype-distinct {shlex.quote(str(params.haplotype_distinct))}"
                )
            if input.highlight:
                parts.append(
                    f"--highlight-file {shlex.quote(str(input.highlight[0]))}"
                )
            shell(" ".join(parts))
