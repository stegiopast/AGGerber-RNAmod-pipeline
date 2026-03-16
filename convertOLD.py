import argparse
import bisect
import logging
from collections import defaultdict

#External dependencies: polars
import polars as pl


logger = logging.getLogger("convert_bedrmod")


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%H:%M:%S",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert transcriptomic BED-like positions to genomic coordinates using a GTF transcript model."
    )
    parser.add_argument("--gtf", required=True, help="Path to GTF annotation file.")
    parser.add_argument("--bed", required=True, help="Path to bedRmod/bed-like transcriptomic file.")
    parser.add_argument(
        "--output",
        default=None,
        help="Optional output TSV path with added columns: enst_id, genomic_pos, map_status.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of ENST IDs to show in unmappable summary table (default: 20).",
    )
    parser.add_argument(
        "--example-n",
        type=int,
        default=10,
        help="Number of unmappable example rows to print (default: 10).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    return parser.parse_args()


def parse_gtf(gtf_path: str) -> pl.DataFrame:
    cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    df = pl.read_csv(
        gtf_path,
        separator="\t",
        has_header=False,
        new_columns=cols,
        comment_prefix="#",
        null_values=".",
        infer_schema_length=1000,
    )
    return df.with_columns(
        [
            pl.col("start").cast(pl.Int64, strict=False),
            pl.col("end").cast(pl.Int64, strict=False),
            pl.col("attributes").str.extract(r'transcript_id\s+"([^"]+)"', 1).alias("transcript_id"),
        ]
    )


def build_transcript_index(gtf_df: pl.DataFrame) -> dict:
    transcript_rows = (
        gtf_df.filter((pl.col("feature") == "transcript") & pl.col("transcript_id").is_not_null())
        .select(["transcript_id", "seqname", "strand"])
        .unique(subset=["transcript_id"])
    )
    transcript_meta = {
        row["transcript_id"]: (row["seqname"], row["strand"])
        for row in transcript_rows.iter_rows(named=True)
    }

    exon_rows = (
        gtf_df.filter((pl.col("feature") == "exon") & pl.col("transcript_id").is_not_null())
        .select(["transcript_id", "start", "end"])
        .iter_rows(named=True)
    )

    exons_by_transcript = defaultdict(list)
    for row in exon_rows:
        if row["start"] is None or row["end"] is None:
            continue
        exon_start_1based = int(row["start"])
        exon_end_1based_inclusive = int(row["end"])
        exon_start_0based = exon_start_1based - 1
        exon_end_0based_exclusive = exon_end_1based_inclusive
        exons_by_transcript[row["transcript_id"]].append((exon_start_0based, exon_end_0based_exclusive))

    transcript_index = {}
    for transcript_id, exon_list in exons_by_transcript.items():
        meta = transcript_meta.get(transcript_id)
        if meta is None:
            continue

        ref_chrom, ref_strand = meta
        exon_list.sort(key=lambda x: x[0], reverse=(ref_strand == "-"))

        starts = []
        ends = []
        cum_lengths = []
        running_total = 0

        for exon_start, exon_end in exon_list:
            exon_len = exon_end - exon_start
            if exon_len <= 0:
                continue
            starts.append(exon_start)
            ends.append(exon_end)
            running_total += exon_len
            cum_lengths.append(running_total)

        if not cum_lengths:
            continue

        transcript_index[transcript_id] = (ref_chrom, ref_strand, starts, ends, cum_lengths)

    return transcript_index


def read_bed_rmod(path: str, has_plain_header: bool = False) -> pl.DataFrame:
    df = pl.read_csv(path, comment_prefix="#", separator="\t", has_header=False)
    if df.width < 2:
        raise ValueError("bedRmod file must contain at least two columns: chrom and chromStart")

    if has_plain_header and df.height > 0:
        df = df.slice(1)

    new_columns = [f"column_{i + 1}" for i in range(df.width)]
    new_columns[0] = "chrom"
    new_columns[1] = "chromStart"
    df.columns = new_columns

    if df["chromStart"].dtype != pl.Int64:
        df = df.with_columns(pl.col("chromStart").cast(pl.Int64, strict=False))

    return df


def read_header_lines(path: str) -> tuple[list[str], str | None]:
    header_lines = []
    plain_header_line = None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
            else:
                candidate = line.rstrip("\n")
                fields = candidate.split("\t")
                if len(fields) >= 2 and fields[0].strip().lower() == "chrom" and fields[1].strip().lower() == "chromstart":
                    plain_header_line = candidate
                break
    return header_lines, plain_header_line


def write_bed_rmod_with_header(
    df: pl.DataFrame, output_path: str, header_lines: list[str], plain_header_line: str | None
) -> None:
    with open(output_path, "w", encoding="utf-8") as handle:
        for hline in header_lines:
            handle.write(hline + "\n")
        if plain_header_line is not None:
            handle.write(plain_header_line + "\n")
        for row in df.iter_rows():
            handle.write("\t".join("" if val is None else str(val) for val in row) + "\n")


def extract_enst(chrom_field: str) -> str | None:
    stringlist = chrom_field.split("|")
    enst_list = [s for s in stringlist if "ENST" in s and "." in s]
    if not enst_list:
        return None
    return enst_list[0]


def convert(transcript_id: str, transcript_pos: int, transcript_index: dict) -> tuple[tuple[str, int, str] | None, str | None]:
    cached = transcript_index.get(transcript_id)
    if cached is None:
        return None, "transcript_not_in_gtf"

    ref_chrom, ref_strand, starts, ends, cum_lengths = cached
    if transcript_pos < 0 or transcript_pos >= int(cum_lengths[-1]):
        return None, "position_out_of_range"

    exon_idx = bisect.bisect_right(cum_lengths, transcript_pos)
    previous_total = 0 if exon_idx == 0 else int(cum_lengths[exon_idx - 1])
    offset = transcript_pos - previous_total

    if ref_strand == "+":
        current_ref_pos = int(starts[exon_idx] + offset) + 1
    else:
        current_ref_pos = int(ends[exon_idx] - 1 - offset) + 1

    return (ref_chrom, current_ref_pos - 1, ref_strand), None


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    logger.info("Reading GTF file: %s", args.gtf)
    gtf_df = parse_gtf(args.gtf)
    logger.info("GTF file read successfully.")

    header_lines, plain_header_line = read_header_lines(args.bed)
    bedrmod = read_bed_rmod(args.bed, has_plain_header=(plain_header_line is not None))
    logger.info("Successfully read bedRmod file with %d entries", bedrmod.shape[0])

    logger.info("Building transcript index")
    transcript_index = build_transcript_index(gtf_df)
    logger.info("Transcript index done")
    logger.info("Indexed transcripts: %d", len(transcript_index))

    gtf_enst_unversioned = {tid.split(".")[0] for tid in transcript_index.keys()}

    total_rows = bedrmod.height
    rows_with_enst = 0
    converted_rows = 0

    unmapped_transcript_not_in_gtf = 0
    unmapped_position_out_of_range = 0
    unmapped_invalid_position = 0
    unmapped_no_enst_in_chrom = 0
    unmapped_version_mismatch = 0
    unmapped_truly_absent = 0

    unmapped_examples = []
    unmapped_by_enst = defaultdict(
        lambda: {"total": 0, "transcript_not_in_gtf": 0, "position_out_of_range": 0, "invalid_position": 0}
    )

    output_chrom = []
    output_chrom_start = []
    output_chrom_end = []
    has_column_3 = "column_3" in bedrmod.columns
    mapped_position_counts = defaultdict(int)

    logger.info("Start converting bedRmod file to genome coordinates")
    iterator = bedrmod.iter_rows(named=True)
    for row in iterator:
        chrom_field = row["chrom"]
        position = row["chromStart"]

        enst_id = extract_enst(chrom_field)
        genomic_pos = None
        status = "no_enst_in_chrom"
        mapped_coords = None

        if enst_id is not None:
            rows_with_enst += 1

            if position is None:
                status = "invalid_position"
                unmapped_invalid_position += 1
                unmapped_by_enst[enst_id]["total"] += 1
                unmapped_by_enst[enst_id]["invalid_position"] += 1
            else:
                converted_tuple, reason = convert(enst_id, int(position), transcript_index=transcript_index)
                if converted_tuple is not None:
                    genomic_chrom, genomic_start_0based, genomic_strand = converted_tuple
                    converted_rows += 1
                    status = "mapped"
                    mapped_coords = (genomic_chrom, genomic_start_0based, genomic_strand)
                    genomic_pos = f"{genomic_chrom}:{genomic_start_0based + 1}:{genomic_strand}"
                    mapped_position_counts[(genomic_chrom, genomic_start_0based, genomic_strand)] += 1
                else:
                    status = reason
                    unmapped_by_enst[enst_id]["total"] += 1

                    if reason == "transcript_not_in_gtf":
                        unmapped_transcript_not_in_gtf += 1
                        unmapped_by_enst[enst_id]["transcript_not_in_gtf"] += 1
                        if enst_id.split(".")[0] in gtf_enst_unversioned:
                            unmapped_version_mismatch += 1
                        else:
                            unmapped_truly_absent += 1
                    elif reason == "position_out_of_range":
                        unmapped_position_out_of_range += 1
                        unmapped_by_enst[enst_id]["position_out_of_range"] += 1

                    if len(unmapped_examples) < args.example_n:
                        unmapped_examples.append((enst_id, int(position), reason))
        else:
            unmapped_no_enst_in_chrom += 1

        if status == "mapped" and mapped_coords is not None:
            mapped_chrom, mapped_start_0based, _ = mapped_coords
            output_chrom.append(mapped_chrom)
            output_chrom_start.append(mapped_start_0based)
            if has_column_3:
                output_chrom_end.append(mapped_start_0based + 1)
        else:
            output_chrom.append(chrom_field)
            output_chrom_start.append(position)
            if has_column_3:
                output_chrom_end.append(row["column_3"])

    if total_rows > 0:
        logger.info(
            "Converted positions (all rows): %d/%d (%.4f)",
            converted_rows,
            total_rows,
            converted_rows / total_rows,
        )

    if rows_with_enst > 0:
        logger.info(
            "Converted positions (rows with ENST): %d/%d (%.4f)",
            converted_rows,
            rows_with_enst,
            converted_rows / rows_with_enst,
        )
    else:
        logger.info("Converted positions (rows with ENST): 0/0 (0.0000)")

    unmapped_total = total_rows - converted_rows
    if unmapped_total > 0:
        logger.warning(
            "Unmappable breakdown: transcript_not_in_gtf=%d, position_out_of_range=%d, no_enst_in_chrom=%d, total=%d",
            unmapped_transcript_not_in_gtf,
            unmapped_position_out_of_range,
            unmapped_no_enst_in_chrom,
            unmapped_total,
        )
        if unmapped_transcript_not_in_gtf > 0:
            logger.warning(
                "Missing transcript details: version_mismatch=%d, truly_absent=%d",
                unmapped_version_mismatch,
                unmapped_truly_absent,
            )
        if unmapped_invalid_position > 0:
            logger.warning("Invalid chromStart (null after cast): %d", unmapped_invalid_position)
        if unmapped_examples:
            logger.warning("First unmappable examples (ENST, chromStart, reason):")
            for ex in unmapped_examples:
                logger.warning("%s", ex)

    if unmapped_by_enst and args.top_n > 0:
        ranked = sorted(unmapped_by_enst.items(), key=lambda x: x[1]["total"], reverse=True)[: args.top_n]
        logger.info("Top %d ENST with most unmappable positions:", len(ranked))
        logger.info("%-20s %10s %12s %14s %12s", "ENST", "total", "not_in_gtf", "out_of_range", "invalid_pos")
        for enst, stats in ranked:
            logger.info(
                "%-20s %10d %12d %14d %12d",
                enst,
                stats["total"],
                stats["transcript_not_in_gtf"],
                stats["position_out_of_range"],
                stats["invalid_position"],
            )

    duplicated_mapped_rows = sum(count - 1 for count in mapped_position_counts.values() if count > 1)
    duplicated_genomic_sites = sum(1 for count in mapped_position_counts.values() if count > 1)
    if converted_rows > 0:
        logger.info(
            "Duplicated genomic positions (mapped rows): %d/%d (%.4f)",
            duplicated_mapped_rows,
            converted_rows,
            duplicated_mapped_rows / converted_rows,
        )
    else:
        logger.info("Duplicated genomic positions (mapped rows): 0/0 (0.0000)")
    logger.info("Genomic sites with duplicates: %d", duplicated_genomic_sites)

    if args.output:
        out_columns = [
            pl.Series(name="chrom", values=output_chrom, strict=False),
            pl.Series(name="chromStart", values=output_chrom_start, strict=False),
        ]
        if has_column_3:
            out_columns.append(pl.Series(name="column_3", values=output_chrom_end, strict=False))

        out_df = bedrmod.with_columns(out_columns)
        ordered_columns = ["chrom", "chromStart"] + [c for c in out_df.columns if c not in {"chrom", "chromStart"}]
        out_df = out_df.select(ordered_columns)

        write_bed_rmod_with_header(out_df, args.output, header_lines, plain_header_line)
        logger.info("Wrote output: %s", args.output)


if __name__ == "__main__":
    main()
