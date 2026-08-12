#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from src import version
from src.utils import build_sample_table, read_sample_fastq_inputs, read_sample_map, stage_fastqs

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MIXCR_LICENSE_FILE = PIPELINE_ROOT / "config" / "mi.license"
DEFAULT_CONTAINERS_FILE = PIPELINE_ROOT / "config" / "containers.json"


def parse_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = value.strip().lower()
    if normalized in {"true", "t", "yes", "y", "1"}:
        return True
    if normalized in {"false", "f", "no", "n", "0"}:
        return False
    raise argparse.ArgumentTypeError("Expected true or false")


def copytree_once(source: Path, target: Path) -> None:
    if not target.exists():
        shutil.copytree(source, target)


def resolve_mixcr_license_file(value: str) -> str:
    if value:
        path = Path(value).expanduser().resolve()
    else:
        path = DEFAULT_MIXCR_LICENSE_FILE
    if path.exists() and path.stat().st_size > 0:
        return str(path)
    return ""


def load_container_manifest() -> dict[str, str]:
    with DEFAULT_CONTAINERS_FILE.open() as handle:
        manifest = json.load(handle)
    required = {"name", "version", "sif", "uri"}
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"Container manifest is missing required keys: {', '.join(missing)}")
    return manifest


def find_container_runtime() -> str:
    for executable in ("singularity", "apptainer"):
        path = shutil.which(executable)
        if path:
            return path
    raise RuntimeError("Neither singularity nor apptainer is available in PATH")


def cache_container(sif_cache: str | Path, dry_run: bool = False) -> Path:
    manifest = load_container_manifest()
    cache_dir = Path(sif_cache).expanduser().resolve()
    sif_path = cache_dir / manifest["sif"]
    if sif_path.exists():
        print(f"Using cached container: {sif_path}")
        return sif_path

    uri = manifest["uri"]
    if dry_run:
        print(f"Would pull {uri} -> {sif_path}")
        return sif_path

    cache_dir.mkdir(parents=True, exist_ok=True)
    runtime = find_container_runtime()
    print(f"Pulling {uri} -> {sif_path}")
    subprocess.run([runtime, "pull", "--force", str(sif_path), uri], check=True)
    return sif_path


def resolve_container_image(args: argparse.Namespace) -> str:
    if args.container_image:
        path = Path(args.container_image).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Container image does not exist: {path}")
        return str(path)
    if args.sif_cache:
        return str(cache_container(args.sif_cache, dry_run=args.dry_run))
    return ""


def resolve_input_fastqs(args: argparse.Namespace) -> list[str]:
    if args.sample_fastq and args.input:
        raise ValueError("Use either --sample-fastq or --input, not both")
    if args.sample_fastq:
        return read_sample_fastq_inputs(args.sample_fastq)
    if args.input:
        return args.input
    raise ValueError("One of --sample-fastq or --input is required")


def resolve_sample_map_path(args: argparse.Namespace) -> str:
    if args.sample_map:
        return args.sample_map
    if args.sample_fastq:
        return args.sample_fastq
    return ""


def init_output(args: argparse.Namespace) -> Path:
    outdir = Path(args.output).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    for name in ["workflow", "resources"]:
        copytree_once(PIPELINE_ROOT / name, outdir / name)
    config_dir = outdir / "config"
    config_dir.mkdir(exist_ok=True)

    input_fastqs = resolve_input_fastqs(args)
    sample_map_path = resolve_sample_map_path(args)
    sample_map = read_sample_map(sample_map_path)
    staged_fastqs = stage_fastqs(input_fastqs, outdir / "input_fastqs", sample_map=sample_map)
    samples = build_sample_table(staged_fastqs, preserve_samples=set(sample_map.values()))
    mixcr_license_file = resolve_mixcr_license_file(args.mixcr_license_file)
    container_image = resolve_container_image(args)
    config = {
        "pipeline": {
            "name": "tcr-seek",
            "version": version,
            "source": str(PIPELINE_ROOT),
        },
        "samples": samples,
        "inputs": {
            "original_fastqs": [str(Path(item).expanduser().absolute()) for item in input_fastqs],
            "staged_fastqs": staged_fastqs,
            "sample_fastq": str(Path(args.sample_fastq).expanduser().absolute()) if args.sample_fastq else "",
            "sample_map": str(Path(sample_map_path).expanduser().absolute()) if sample_map_path else "",
            "sample_info": str(Path(args.sample_info).expanduser().absolute()) if args.sample_info else "",
        },
        "paths": {
            "output": str(outdir),
            "refdir": str(Path(args.refdir).resolve()),
            "metadata": str(Path(args.metadata).resolve()) if args.metadata else "",
            "sample_info": str(Path(args.sample_info).resolve()) if args.sample_info else "",
            "mixcr_license_file": mixcr_license_file,
            "container_image": container_image,
            "container_manifest": str(DEFAULT_CONTAINERS_FILE),
            "sif_cache": str(Path(args.sif_cache).expanduser().resolve()) if args.sif_cache else "",
        },
        "params": {
            "threads": int(args.threads),
            "species": args.species,
            "receptor": args.receptor,
            "format": "airr",
            "clonotype_definition": args.clonotype_definition,
            "run_rarefaction": bool(args.run_rarefaction),
            "cleanup_intermediates": bool(args.cleanup_intermediates),
            "qc_variables": args.qc_variables,
        },
    }
    with (outdir / "config.json").open("w") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
    return outdir


def base_snakemake_cmd(outdir: Path, printshellcmds: bool = True) -> list[str]:
    cmd = [
        "snakemake",
        "-s",
        str(outdir / "workflow" / "Snakefile"),
        "-d",
        str(outdir),
        "--configfile",
        str(outdir / "config.json"),
        "--rerun-incomplete",
        "--keep-going",
    ]
    if printshellcmds:
        cmd.append("--printshellcmds")
    return cmd


def run_local(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> None:
    cmd = base_snakemake_cmd(outdir) + [
        "--cores",
        str(args.cores),
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    subprocess.run(cmd, check=True, env=env)


def slurm_dirs(outdir: Path) -> tuple[Path, Path]:
    logdir = outdir / "logfiles"
    slurm_dir = logdir / "slurm"
    logdir.mkdir(exist_ok=True)
    slurm_dir.mkdir(exist_ok=True)
    return logdir, slurm_dir


def slurm_snakemake_cmd(args: argparse.Namespace, outdir: Path, slurm_dir: Path, printshellcmds: bool = True) -> list[str]:
    cluster_cmd = " ".join(
        [
            "sbatch",
            "--parsable",
            "--partition={resources.partition}",
            "--time={resources.walltime}",
            "--mem={resources.mem_mb}",
            "--cpus-per-task={threads}",
            "--job-name=tcr-seek.{rule}",
            f"--output={shlex.quote(str(slurm_dir))}/slurm-%j_{{rule}}.out",
            f"--error={shlex.quote(str(slurm_dir))}/slurm-%j_{{rule}}.out",
        ]
    )
    return base_snakemake_cmd(outdir, printshellcmds=printshellcmds) + [
        "--jobs",
        str(args.jobs),
        "--latency-wait",
        str(args.latency_wait),
        "--restart-times",
        str(args.restart_times),
        "--local-cores",
        str(args.local_cores),
        "--use-singularity",
        "--singularity-args",
        f"--bind {args.singularity_bind}",
        "--cluster",
        cluster_cmd,
    ]


def write_slurm_controller_script(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> Path:
    logdir, slurm_dir = slurm_dirs(outdir)
    cmd = slurm_snakemake_cmd(args, outdir, slurm_dir, printshellcmds=True)
    script = outdir / "run_tcr_seek_controller.sh"
    exports = []
    if env.get("MI_LICENSE_FILE"):
        exports.append(f"export MI_LICENSE_FILE={shlex.quote(env['MI_LICENSE_FILE'])}")
    exports_text = "\n".join(exports)
    script.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "if command -v module >/dev/null 2>&1; then\n"
        "  module load snakemake singularity\n"
        "fi\n\n"
        f"cd {shlex.quote(str(outdir))}\n"
        f"{exports_text}\n\n"
        f"{shlex.join(cmd)}\n"
    )
    script.chmod(0o755)
    return script


def run_slurm_dry_run(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> None:
    logdir, slurm_dir = slurm_dirs(outdir)
    cmd = slurm_snakemake_cmd(args, outdir, slurm_dir, printshellcmds=not args.quiet)
    cmd.append("--dry-run")
    log_file = logdir / "dry-run.log"
    result = subprocess.run(cmd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log_file.write_text(result.stdout)
    if not args.quiet:
        print(result.stdout, end="")
    if result.returncode != 0:
        print(f"Dry-run failed. See: {log_file}", file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout)
    print(f"Dry-run complete. Full Snakemake plan written to: {log_file}")


def submit_slurm_controller(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> None:
    logdir, _ = slurm_dirs(outdir)
    script = write_slurm_controller_script(args, outdir, env)
    stdout = logdir / "controller-%j.out"
    submit_cmd = [
        "sbatch",
        "--parsable",
        "--partition",
        args.controller_partition,
        "--time",
        args.controller_time,
        "--mem",
        str(args.controller_mem_mb),
        "--cpus-per-task",
        str(args.local_cores),
        "--job-name",
        args.job_name,
        "--output",
        str(stdout),
        "--error",
        str(stdout),
        str(script),
    ]
    result = subprocess.run(submit_cmd, check=True, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    job_id = result.stdout.strip()
    print(f"Submitted batch job {job_id}")
    print(f"Controller log: {stdout}")
    print(f"Worker logs: {outdir / 'logfiles' / 'slurm'}")


def run_slurm_direct(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> None:
    _, slurm_dir = slurm_dirs(outdir)
    cmd = slurm_snakemake_cmd(args, outdir, slurm_dir, printshellcmds=True)
    subprocess.run(cmd, check=True, env=env)


def run_slurm(args: argparse.Namespace, outdir: Path, env: dict[str, str]) -> None:
    if not (args.container_image or args.sif_cache):
        raise ValueError("Slurm mode requires --container-image or --sif-cache")
    if args.dry_run:
        run_slurm_dry_run(args, outdir, env)
    elif args.no_submit_controller:
        run_slurm_direct(args, outdir, env)
    else:
        submit_slurm_controller(args, outdir, env)


def run_pipeline(args: argparse.Namespace) -> None:
    outdir = init_output(args)
    env = os.environ.copy()
    mixcr_license_file = resolve_mixcr_license_file(args.mixcr_license_file)
    if mixcr_license_file:
        env["MI_LICENSE_FILE"] = mixcr_license_file
    if args.mode == "local":
        run_local(args, outdir, env)
    elif args.mode == "slurm":
        run_slurm(args, outdir, env)
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")


def cache_command(args: argparse.Namespace) -> None:
    cache_container(args.sif_cache, dry_run=args.dry_run)


def unlock(args: argparse.Namespace) -> None:
    outdir = Path(args.output).resolve()
    cmd = [
        "snakemake",
        "--unlock",
        "-s",
        str(outdir / "workflow" / "Snakefile"),
        "-d",
        str(outdir),
        "--configfile",
        str(outdir / "config.json"),
        "--cores",
        "1",
    ]
    subprocess.run(cmd, check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="tcr-seek", description="Bulk TCR-seq processing and immunarch import pipeline")
    parser.add_argument("--version", action="version", version=f"tcr-seek {version}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run pRESTO, Change-O/IgBLAST, immunarch object creation, and diversity metrics")
    run.add_argument("--input", nargs="+", help="Paired FASTQ.gz files. Names are normalized to cell-seek style sample.R1.fastq.gz/sample.R2.fastq.gz symlinks.")
    run.add_argument("--sample-fastq", default="", help="TSV with fastq and sample columns. The fastq column supplies inputs, and the same file is used as --sample-map unless --sample-map is provided.")
    run.add_argument("--output", required=True, help="Pipeline output directory")
    run.add_argument("--refdir", required=True, help="Directory with IS_Human_R1_Primers.fasta, IS_Human_R2_Primers.fasta, IS_Human_C-Region.fasta, and Immune_Ref.fasta")
    run.add_argument("--metadata", default="", help="Optional tab-delimited metadata with a Sample column")
    run.add_argument("--sample-info", default="", help="Optional tab-delimited sample information file. Must contain Sample, SampleID, sample, or sample_id; columns are merged into immdata metadata.")
    run.add_argument("--sample-map", default="", help="Optional TSV with fastq and sample columns. Listed FASTQs use the provided sample name; unlisted FASTQs use automatic cell-seek-style parsing.")
    run.add_argument("--mixcr-license-file", default="", help="Optional path to a MiXCR mi.license file. Defaults to tcr-seek/config/mi.license when that file is non-empty.")
    run.add_argument("--species", default="human", help="Species for changeo-igblast, default: human")
    run.add_argument("--receptor", default="tr", help="Receptor type for changeo-igblast, default: tr")
    run.add_argument("--clonotype-definition", default="strict", help="Clonotype definition preset or comma-separated immunarch columns for diversity metrics. Presets: gene, nt, aa, strict. Default: strict")
    run.add_argument("--run-rarefaction", type=parse_bool, default=True, metavar="{true,false}", help="Whether to calculate rarefaction metrics and curves for the QC report. Default: true")
    run.add_argument("--cleanup-intermediates", type=parse_bool, default=False, metavar="{true,false}", help="After successful final outputs, gzip retained Presto handoff FASTQs and remove reproducible intermediates. Default: false")
    run.add_argument("--qc-variables", default="", help="Comma-separated sample metadata variables to use as QC report grouping tabs. Default: all non-Sample sample-info/metadata columns.")
    run.add_argument("--threads", type=int, default=2, help="Threads per sample for pRESTO steps")
    run.add_argument("--cores", type=int, default=2, help="Total Snakemake local cores")
    run.add_argument("--mode", choices=["local", "slurm"], default="local", help="Execution mode. local runs in the current process; slurm submits rule jobs to the scheduler.")
    run.add_argument("--container-image", default="", help="Explicit SIF path used to run Slurm worker jobs. Overrides --sif-cache.")
    run.add_argument("--sif-cache", default="", help="Directory containing or receiving the default tcr-seek SIF. Used when --container-image is not provided.")
    run.add_argument("--jobs", type=int, default=100, help="Maximum number of concurrent Slurm jobs submitted by Snakemake")
    run.add_argument("--local-cores", type=int, default=2, help="Cores reserved for the Snakemake controller in Slurm mode")
    run.add_argument("--job-name", default="tcr-seek.controller", help="Slurm job name for the tcr-seek controller job")
    run.add_argument("--controller-partition", default="norm", help="Slurm partition for the tcr-seek controller job")
    run.add_argument("--controller-time", default="1-00:00:00", help="Walltime for the tcr-seek controller job")
    run.add_argument("--controller-mem-mb", type=int, default=8000, help="Memory in MB for the tcr-seek controller job")
    run.add_argument("--no-submit-controller", action="store_true", help="Run the Slurm Snakemake controller in the foreground instead of submitting a controller job")
    run.add_argument("--quiet", action="store_true", help="Do not print captured dry-run output to the terminal; still writes it to logfiles/dry-run.log")
    run.add_argument("--latency-wait", type=int, default=120, help="Seconds Snakemake waits for files to appear on shared filesystems")
    run.add_argument("--restart-times", type=int, default=1, help="Number of times to retry failed Slurm jobs")
    run.add_argument("--singularity-bind", default="/data:/data,/vf:/vf,/tmp:/tmp", help="Comma-separated Singularity bind paths for Slurm worker jobs")
    run.add_argument("--dry-run", action="store_true", help="Create output/config and dry-run Snakemake")
    run.set_defaults(func=run_pipeline)

    cache_parser = sub.add_parser("cache", help="Pull the default tcr-seek container into a shareable SIF cache")
    cache_parser.add_argument("--sif-cache", required=True, help="Directory where the default tcr-seek SIF will be stored")
    cache_parser.add_argument("--dry-run", action="store_true", help="Show the container pull that would be performed")
    cache_parser.set_defaults(func=cache_command)

    unlock_parser = sub.add_parser("unlock", help="Unlock an existing output directory")
    unlock_parser.add_argument("--output", required=True)
    unlock_parser.set_defaults(func=unlock)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
