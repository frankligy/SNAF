#!/usr/bin/env python3
import argparse
import os
import subprocess
from datetime import datetime
from typing import Tuple

CONTAINER_IMAGE = "/data/salomonis2/LabFiles/Frank-Li/neoantigen/alphafold/alphafold.sif"
ROOT_MOUNT_DIRECTORY = "/mnt"


def main():
    args = parse_arguments()

    # Path to the Uniref90 database for use by JackHMMER.
    uniref90_database_path = os.path.join(args.data_dir, "uniref90", "uniref90.fasta")

    # Path to the Uniprot database for use by JackHMMER.
    uniprot_database_path = os.path.join(args.data_dir, "uniprot", "uniprot.fasta")

    # Path to the MGnify database for use by JackHMMER.
    mgnify_database_path = os.path.join(
        args.data_dir, "mgnify", "mgy_clusters_2018_12.fa"
    )

    # Path to the BFD database for use by HHblits.
    bfd_database_path = os.path.join(
        args.data_dir, "bfd", "bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    )

    # Path to the Small BFD database for use by JackHMMER.
    small_bfd_database_path = os.path.join(
        args.data_dir, "small_bfd", "bfd-first_non_consensus_sequences.fasta"
    )

    # Path to the Uniclust30 database for use by HHblits.
    uniclust30_database_path = os.path.join(
        args.data_dir, "uniclust30", "uniclust30_2018_08", "uniclust30_2018_08"
    )

    # Path to the PDB70 database for use by HHsearch.
    pdb70_database_path = os.path.join(args.data_dir, "pdb70", "pdb70")

    # Path to the PDB seqres database for use by hmmsearch.
    pdb_seqres_database_path = os.path.join(
        args.data_dir, "pdb_seqres", "pdb_seqres.txt"
    )

    # Path to a directory with template mmCIF structures, each named <pdb_id>.cif')
    template_mmcif_dir = os.path.join(args.data_dir, "pdb_mmcif", "mmcif_files")

    # Path to a file mapping obsolete PDB IDs to their replacements.
    obsolete_pdbs_path = os.path.join(args.data_dir, "pdb_mmcif", "obsolete.dat")

    mounts = []
    command_args = []

    # Mount each fasta path as a unique target directory
    target_fasta_paths = []
    for i, fasta_path in enumerate(args.fasta_paths):
        mount, target_path = _generate_mount(f"fasta_path_{i}", fasta_path)
        mounts.append(mount)
        target_fasta_paths.append(target_path)
    command_args.append(f"--fasta_paths={','.join(target_fasta_paths)}")

    # Mount database and output directories
    database_paths = [
        ("uniref90_database_path", uniref90_database_path),
        ("mgnify_database_path", mgnify_database_path),
        ("data_dir", args.data_dir),
        ("template_mmcif_dir", template_mmcif_dir),
        ("obsolete_pdbs_path", obsolete_pdbs_path),
    ]
    if args.model_preset == "multimer":
        database_paths.append(("uniprot_database_path", uniprot_database_path))
        database_paths.append(("pdb_seqres_database_path", pdb_seqres_database_path))
    else:
        database_paths.append(("pdb70_database_path", pdb70_database_path))

    if args.db_preset == "reduced_dbs":
        database_paths.append(("small_bfd_database_path", small_bfd_database_path))
    else:
        database_paths.extend(
            [
                ("uniclust30_database_path", uniclust30_database_path),
                ("bfd_database_path", bfd_database_path),
            ]
        )

    for name, path in database_paths:
        if path:
            mount, target_path = _generate_mount(name, path)
            mounts.append(mount)
            command_args.append(f"--{name}={target_path}")

    output_mount, output_target_path = _generate_mount(
        "output", args.output_dir, read_only=False
    )
    mounts.append(output_mount)

    # Set general options for the alphafold script
    command_args.extend(
        [
            f"--output_dir={output_target_path}",
            f"--max_template_date={args.max_template_date}",
            f"--db_preset={args.db_preset}",
            f"--model_preset={args.model_preset}",
            f"--benchmark={args.benchmark}",
            f"--use_precomputed_msas={args.use_precomputed_msas}",
            "--logtostderr",
        ]
    )

    # Set environment variables for the container
    env = {
        "NVIDIA_VISIBLE_DEVICES": args.gpu_devices,
        # The following flags allow us to make predictions on proteins that
        # would typically be too long to fit into GPU memory.
        "TF_FORCE_UNIFIED_MEMORY": "1",
        "XLA_PYTHON_CLIENT_MEM_FRACTION": "4.0",
        "MAX_CPUS": args.cpus,
    }

    # Generate the final command to execute
    command = [
        "singularity",
        "exec",
        "--nv" if args.use_gpu else "",
        "--bind",
        ",".join(mounts),
        *[f'--env="{k}={v}"' for k, v in env.items()],
        CONTAINER_IMAGE,
        "/app/run_alphafold.sh",
        *command_args,
    ]

    print("Executing: " + " ".join(command))

    p = subprocess.run(command)
    p.check_returncode()


def _generate_mount(mount_name: str, path: str, read_only=True) -> Tuple[str, str]:
    """
    Generate a mount line for a singularity container.
    :param mount_name: The name of the mount point.
    :param path: The path to mount.
    :return: A tuple of the mount line and the path to mount.
    """
    path = os.path.abspath(path)
    source_path = os.path.dirname(path)
    target_path = os.path.join(ROOT_MOUNT_DIRECTORY, mount_name)
    opts = "ro" if read_only else "rw"

    mount_cmd = f"{source_path}:{target_path}:{opts}"
    return mount_cmd, os.path.join(target_path, os.path.basename(path))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Singularity launch script for Alphafold v2.1.0"
    )

    parser.add_argument(
        "--fasta-paths",
        "-f",
        required=True,
        nargs="+",
        help="Paths to FASTA files, each containing one sequence. "
        "All FASTA paths must have a unique basename as the basename "
        "is used to name the output directories for each prediction.",
    )
    parser.add_argument(
        "--is-prokaryote-list",
        required=False,
        nargs="+",
        help="Optional for multimer system, "
        "not used by the single chain system. "
        "This list should contain a boolean for each fasta "
        "specifying true where the target complex is from a "
        "prokaryote, and false where it is not, or where the "
        "origin is unknown. These values determine the pairing "
        "method for the MSA.",
    )
    parser.add_argument(
        "--max-template-date",
        "-t",
        default=datetime.today().strftime("%Y-%m-%d"),
        help="Maximum template release date to consider "
        "(ISO-8601 format - i.e. YYYY-MM-DD). "
        "Important if folding historical test sets.",
    )
    parser.add_argument(
        "--db-preset",
        choices=["reduced_dbs", "full_dbs"],
        default="full_dbs",
        help="Choose preset model configuration - no ensembling with "
        "uniref90 + bfd + uniclust30 (full_dbs), or "
        "8 model ensemblings with uniref90 + bfd + uniclust30 (casp14).",
    )
    parser.add_argument(
        "--model-preset",
        choices=["monomer", "monomer_casp14", "monomer_ptm", "multimer"],
        default="monomer",
        help="Choose preset model configuration - the monomer model, the monomer model "
        "with extra ensembling, monomer model with pTM head, or multimer model",
    )
    parser.add_argument(
        "--benchmark",
        "-b",
        default=False,
        action="store_true",
        help="Run multiple JAX model evaluations to obtain a timing "
        "that excludes the compilation time, which should be more indicative "
        "of the time required for inferencing many proteins.",
    )
    parser.add_argument(
        "--use-precomputed-msas",
        default=False,
        action="store_true",
        help="Whether to read MSAs that have been written to disk. WARNING: This will "
        "not check if the sequence, database or configuration have changed.",
    )
    parser.add_argument(
        "--data-dir",
        "-d",
        default="/database/alphafold/databases",
        help="Path to directory with supporting data: AlphaFold parameters and genetic "
        "and template databases. Set to the target of download_all_databases.sh.",
    )
    parser.add_argument(
        "--docker-image", default=CONTAINER_IMAGE, help="Alphafold docker image."
    )
    parser.add_argument(
        "--output-dir", "-o", default="results/", help="Output directory for results."
    )
    parser.add_argument(
        "--use-gpu",
        default=True,
        action="store_true",
        help="Enable NVIDIA runtime to run with GPUs.",
    )
    parser.add_argument(
        "--gpu-devices",
        default="all",
        help="Comma separated list of devices to pass to NVIDIA_VISIBLE_DEVICES.",
    )
    parser.add_argument(
        "--cpus", "-c", type=int, default=8, help="Number of CPUs to use."
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
