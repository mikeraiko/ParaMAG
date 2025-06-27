# cli_wrapper.py
import argparse
import os
import subprocess
import sys
import yaml


def write_config(reads1, reads2, adapter_file, config_path):
    config = {
        'reads1': reads1,
        'reads2': reads2,
        'adapter_file': adapter_file or 'adapters.fa'
    }
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    print(f"[INFO] Generated config file at {config_path}")


def validate_args(args):
    if not os.path.exists(args.reads1):
        sys.exit(f"Error: Input file {args.reads1} not found.")
    if not os.path.exists(args.reads2):
        sys.exit(f"Error: Input file {args.reads2} not found.")
    if args.adapter_file and not os.path.exists(args.adapter_file):
        sys.exit(f"Error: Adapter file {args.adapter_file} not found.")


def run_pipeline(args):
    config_path = args.config or "autogen_config.yaml"
    write_config(args.reads1, args.reads2, args.adapter_file, config_path)

    snakemake_cmd = [
        "snakemake",
        "--cores", str(args.threads),
        "--use-conda",
        "--snakefile", "snakefile",
        "--configfile", config_path,
        "--rerun-incomplete",
        "--printshellcmds",
        "--latency-wait", "30"
    ]
    print("\n[INFO] Running pipeline...")
    subprocess.run(snakemake_cmd, check=True)
    print("\n[INFO] Pipeline finished successfully.")


def run_example():
    print("[INFO] Downloading and running example dataset...")
  #  os.mkdir("example_run")
  #  os.chdir("example_run")
    subprocess.run(["python", "paramag.py", "--reads1", "example_data/sample_R1.fastq",
                    "--reads2", "example_data/sample_R2.fastq",
                    "--adapter-file", "adapters.fa",
                    "--threads", "4"])




def main():
    parser = argparse.ArgumentParser(description="ParaMAG pipeline wrapper")
    parser.add_argument("--reads1", help="Path to R1 FASTQ file")
    parser.add_argument("--reads2", help="Path to R2 FASTQ file")
    parser.add_argument("--adapter-file", help="Path to adapter file (default: adapters.fa)", default="adapters.fa")
    parser.add_argument("--config", help="Output config file (optional)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--example", action="store_true", help="Run pipeline on included example dataset")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if args.example:
        run_example()
    else:
        validate_args(args)
        run_pipeline(args)


if __name__ == "__main__":
    main()
