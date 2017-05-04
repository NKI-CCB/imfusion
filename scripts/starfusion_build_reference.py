"""Script for building StarFusion indices using FusionFilter."""

import argparse
import gzip
import logging
from pathlib import Path
import shutil
import subprocess

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function."""

    args = parse_args()

    # Create output and temp dirs.
    tmp_dir = args.output_dir / '_tmp'
    tmp_dir.mkdir(parents=True, exist_ok=False)

    # Filter the patch chromosomes from the gtf, as these are
    # likely not present in the fasta.
    logging.info('- Generating cDNA sequences')

    gtf_path = tmp_dir / 'ref.gtf'
    with gtf_path.open('wb') as file_:
        check_call(['grep', '-v', r'^\(MG\|JH\|GL\)', str(args.gtf)],
                   stdout=file_)

    # Create cDNA_seqs file.
    cdna_path = tmp_dir / 'cDNA_seqs.fa'
    with cdna_path.open('wb') as file_:
        script_path = args.ff_path / 'util' / 'gtf_file_to_cDNA_seqs.pl'
        check_call(['perl', str(script_path), str(gtf_path), str(args.fasta)],
                   stdout=file_)

    # Build masked cDNA_seqs file using RepeatMasker.
    # Note: requires library to be installed from http://www.girinst.org.
    logging.info('- Masking repeats')

    masked_path = cdna_path.with_suffix('.fa.masked')
    check_call([str(args.rm_path / 'RepeatMasker'), '-pa', str(args.threads), '-s',
                '-species', 'mouse', '-xsmall', str(cdna_path)])

    # Create blastpairs.
    logging.info('- Creating blast pairs')

    check_call(['makeblastdb', '-in', str(masked_path), '-dbtype', 'nucl'])

    pair_path = tmp_dir / 'blast_pairs.outfmt6'
    with pair_path.open('wb') as file_:
        check_call(['blastn',
                    '-query', str(cdna_path),
                    '-db', str(masked_path),
                    '-max_target_seqs', '10000',
                    '-outfmt', '6',
                    '-evalue', '1e-3',
                    '-lcase_masking',
                    '-num_threads', str(args.threads),
                    '-word_size', '11'],
                   stdout=file_)

    pair_gz_path = pair_path.with_suffix('.gene_syms.outfmt6.gz')
    with gzip.open(str(pair_gz_path), 'wb') as file_:
        script_path = (args.ff_path / 'util' /
                       'blast_outfmt6_replace_trans_id_w_gene_symbol.pl')
        check_call(['perl', str(script_path), str(cdna_path), str(pair_path)],
                   stdout=file_)

    # Prepare library.
    logging.info('- Preparing library')
    script_path = args.ff_path / 'util' / 'prep_genome_lib.pl'
    check_call(['perl', str(script_path),
                '--genome_fa', str(args.fasta),
                '--gtf', str(gtf_path),
                '--blast_pairs', str(pair_gz_path),
                '--cdna_fa', str(cdna_path),
                '--CPU', str(args.threads),
                '--max_readlength', str(args.read_length),
                '--output_dir', str(args.output_dir)])

    # shutil.rmtree(str(args.tmp_dir))

def parse_args():
    """Parses command line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta', required=True, type=Path)
    parser.add_argument('--gtf', required=True, type=Path)
    parser.add_argument('--output_dir', required=True, type=Path)
    
    parser.add_argument('--ff_path', required=False, default='', type=Path)
    parser.add_argument('--rm_path', required=False, default='', type=Path)

    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--read_length', type=int, default=100)

    return parser.parse_args()


def check_call(args, verbose=True, **kwargs):
    """Wrapper function for check_call."""

    if verbose:
        logging.info('Running command - "%s"', ' '.join(args))

    subprocess.check_call(args, **kwargs)


if __name__ == '__main__':
    main()
