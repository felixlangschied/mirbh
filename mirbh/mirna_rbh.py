import glob
import subprocess as sp
import os
import sys
import argparse
import multiprocessing as mp


class Sequence(object):
    def __init__(self, name, chromosome, start, end, strand, seq):
        # chromosome
        if chromosome.startswith('chr'):
            self.chromosome = chromosome.split('chr')[1]
        else:
            self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand
        # miRNA identifier
        self.name = name
        # nucleotide sequence, always as DNA not RNA
        self.seq = seq.replace('U', 'T')
        self.seq = seq.replace('-', '')


def check_blastdb(db_path):
    file_extensions = ['nhr', 'nin', 'nsq']
    for fe in file_extensions:
        files = glob.glob(f'{db_path}*{fe}')
        if not files:
            # At least one of the BLAST db files is not existent
            return False
    return True


def make_blastndb(inpath, outpath):
    db_command = 'makeblastdb -in {} -out {} -dbtype nucl'.format(inpath, outpath)
    sp.call(db_command, shell=True)


def run_blast(seq, db, threads, dust_f='no', evalue=10):
    blast_cmd = (
        'blastn -task blastn -num_threads "{}" -dust "{}" -db {} -evalue {} '
        '-outfmt "6 saccver sstart send sstrand length sseq"'.format(threads, dust_f, db, evalue)
    )
    blast_call = sp.Popen(
        blast_cmd, shell=True, stdin=sp.PIPE,
        stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    seq = seq.replace('-', '')
    # run BLAST
    res, err = blast_call.communicate(seq)
    if err:
        print(f'ERROR: {err}')
        sys.exit()
    # parse BLASTn results
    for result in res.split('\n'):
        # last line of results is empty, therefore we need to skip empty rows
        if result:
            return result
    print('No BLAST hit above length cutoff')
    # best_hit = res.split('\n')[0]
    # return best_hit


def loc_check(mirna, hit):
    # compare re-BLAST hit with reference miRNA
    if (
            mirna.chromosome == hit.chromosome
            and mirna.strand == hit.strand
    ):
        # first within second
        if (
                (hit.start <= mirna.start <= hit.end)
                or (hit.start <= mirna.end <= hit.end)
        ):
            return True
        # second within first
        elif (
                (mirna.start <= hit.start <= mirna.end)
                or (mirna.start <= hit.end <= mirna.end)
        ):
            return True
        # No overlap
        else:
            return False
    else:
        return False


def main():
    ##########################################################################
    # Command line arguments
    ##########################################################################
    parser = argparse.ArgumentParser(
        description='Perform reciprocal best BLAST hit search for reference miRNAs'
    )
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    # covariance models folder
    required.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str, required=True,
        help='Path to a tab seperated file of reference miRNAs'
    )
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str, required=True,
        help='Path to the output directory'
    )
    # query genome
    required.add_argument(
        '-q', '--query', metavar='<.fa>', type=str, required=True,
        help='Path to query genome in FASTA format'
    )
    # reference genome
    required.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str, required=True,
        help='Path to reference genome in FASTA format'
    )
    ##########################################################################
    # Optional Arguments
    ##########################################################################
    # query_name
    optional.add_argument(
        '--queryname', metavar='str', type=str, nargs='?', const='', default='',
        help=(
            'Name the output file'
        )
    )
    # cpu, use maximum number of available cpus unless specified otherwise
    optional.add_argument(
        '--cpu', metavar='int', type=int,
        help='Number of CPU cores to use (Default: all available)', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # length filter to prevent short hits
    optional.add_argument(
        '--minlength', metavar='float', type=float,
        help='Reciprocal hit in the query species must have at '
             'least the length of this value times the length of the refernce pre-miRNA (Default: 0.5)',
        nargs='?', const=0.5, default=0.5
    )
    optional.add_argument(
        '--evalue', metavar='float', type=float,
        help='Minimum Expect value for BLAST hits (Default: 10)',
        nargs='?', const=10, default=10
    )
    # use dust filter?
    parser.add_argument(
        '--dust', metavar='yes/no', type=str,
        help='Use BLASTn dust filter during re-BLAST.',
        nargs='?',
        const='no', default='no'
    )
    ##########################################################################
    # Parse Arguments
    ##########################################################################
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    ref_path = args.reference
    n_path = args.ncrna
    q_path = args.query
    out_dir = args.output
    query_name = args.queryname
    cpu = args.cpu
    dust = args.dust
    length_ratio = args.minlength
    exp = args.evalue

    ##########################################################################
    # Input checks
    ##########################################################################
    if os.path.isfile(ref_path):
        if not check_blastdb(ref_path):
            print('# Making reference database')
            make_blastndb(ref_path, ref_path)
    else:
        print('ERROR: Reference genome not found')
        sys.exit()
    if os.path.isfile(q_path):
        if not check_blastdb(q_path):
            print('# Making query database')
            make_blastndb(q_path, q_path)
    else:
        print('ERROR: Query genome not found')
        sys.exit()
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if cpu > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    if not query_name:
        query_name = q_path.split('/')[-1].split('.')[0]

    ##########################################################################
    # Main algorithm
    ##########################################################################

    print('# Reading miRNA information')
    mirna_dict = {}
    with open(n_path, 'r') as fh:
        for line in fh:
            name, chromosome, start, end, strand, pre, mature = line.strip().split('\t')
            mirna_dict[name] = Sequence(name, chromosome, start, end, strand, pre)
    print('# Done')

    outfile = f'{out_dir}/{query_name}_orthologs.fa'
    with open(outfile, 'w') as of:
        # determine distribution of miRNA length
        for mirid in mirna_dict:
            print(f'### {mirid}')
            mirna = mirna_dict[mirid]
            mirna_length = mirna.end - mirna.start
            length_cutoff = mirna_length * length_ratio
            # BLAST in query genome
            preseq = mirna_dict[mirid].seq
            best_hit = run_blast(preseq, q_path, cpu, dust, exp)
            if best_hit:
                hit_seq = best_hit.split()[-1]
            else:
                print(f'{mirid} not found in query genome')
                continue

            # re-BLAST
            best_rb = run_blast(hit_seq, ref_path, cpu, dust, exp)
            if best_rb:
                hit_chrom, hit_start, hit_end, hit_strand, hit_length, hit_seq = best_rb.split()
                if int(hit_length) < length_cutoff:
                    print(f'REJECTED: Best re-BLAST hit of {mirid} below length threshold')
                    continue
                if hit_strand == 'plus':
                    hit_strand = '+'
                elif hit_strand == 'minus':
                    hit_strand = '-'
                    # switch positions so that smaller location is always start
                    tmp = hit_start
                    hit_start = hit_end
                    hit_end = tmp
            else:
                print(f're-BLAST of {mirid} yielded no results')
                continue

            hit = Sequence(mirid, hit_chrom, hit_start, hit_end, hit_strand, hit_seq)

            # check location of hit
            print(f'mirna: {mirna.start} - {mirna.end}')
            print(f'hit: {hit.start} - {hit.end}')
            if loc_check(mirna, hit):
                print('CONFIRMED')
                outstr = f'>{query_name}|{mirid}|{hit.chromosome}|{hit.start}|{hit.end}|{hit.strand}\n{hit.seq}\n'
                of.write(outstr)
            else:
                print('REJECTED: Position of best re-BLAST hit does not overlap with reference')


if __name__ == "__main__":
    main()
