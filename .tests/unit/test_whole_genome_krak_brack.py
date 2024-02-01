import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_whole_genome_krak_brack():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/whole_genome_krak_brack/data")
        expected_path = PurePosixPath(".tests/unit/whole_genome_krak_brack/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("output/assemblySpeciesID/wgkb/sample1/sample1.kraken_taxonomy.txt output/assemblySpeciesID/wgkb/sample1/sample1.kraken_output.txt output/assemblySpeciesID/wgkb/sample1/sample1.bracken_output.txt", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "output/assemblySpeciesID/wgkb/sample1/sample1.kraken_taxonomy.txt output/assemblySpeciesID/wgkb/sample1/sample1.kraken_output.txt output/assemblySpeciesID/wgkb/sample1/sample1.bracken_output.txt",
            "-f", 
            "-j1",
            "--target-files-omit-workdir-adjustment",
    
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
