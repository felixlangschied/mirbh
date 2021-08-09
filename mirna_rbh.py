import glob
import subprocess as sp
import os
import sys


ref_path = '/share/project2/felix/ncOrtho/mirgenedb/05cutoff/Homo_sapiens/data/Homo_sapiens.fa'
n_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/test.tsv'
# n_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/mirgenedb.tsv'
# q_path = '/share/project2/felix/ncOrtho/mirgenedb/05cutoff/Homo_sapiens/data/Homo_sapiens.fa'
q_path = '/share/project2/felix/ncOrtho/mirgenedb/05cutoff/Mus_musculus/data/Mus_musculus.fa'
out_dir = '/share/project2/felix/ncOrtho/rbh_test'
query_name = 'Mus_musculus'
cpu = '4'
dust = 'no'
length_ratio = 0.7


# Central class of microRNA objects
class Sequence(object):
    def __init__(self, name, chromosome, start, end, strand, seq):
        # chromosome that the miRNA is located on
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
        # nucleotide sequence of the pre-miRNA
        self.seq = seq.replace('U', 'T')


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


def run_blast(seq, db, len_c, threads, dust_f='no'):
    blast_cmd = (
        f'blastn -task blastn -num_threads "{threads}" -dust "{dust_f}" -db {db} -outfmt "6 saccver sstart send sstrand length sseq"'
    )
    blast_call = sp.Popen(
        blast_cmd, shell=True, stdin=sp.PIPE,
        stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    # remove hyphens from sequence
    seq = seq.replace('-', '')
    res, err = blast_call.communicate(seq)
    if err:
        print(f'ERROR: {err}')
        sys.exit()
    # parse BLASTn results
    for result in res.split('\n'):
        # last line of results is empty, therefore skip
        if result:
            length = int(result.split()[-2])
            if length >= len_c:
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


def main():
    # input checks
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

    print('# Reading miRNA information')
    mirna_dict = {}
    with open(n_path, 'r') as fh:
        for line in fh:
            name, chromosome, start, end, strand, pre, mature = line.strip().split('\t')
            mirna_dict[name] = Sequence(name, chromosome, start, end, strand, pre)
    print('# Done')

    outfile = f'{out_dir}/{query_name}.fa'
    with open(outfile, 'w') as of:
        for mirid in mirna_dict:
            print(f'### {mirid}')
            mirna = mirna_dict[mirid]
            mirna_length = mirna.end - mirna.start
            length_cutoff = mirna_length * length_ratio
            # BLAST in query genome
            preseq = mirna_dict[mirid].seq
            best_hit = run_blast(preseq, q_path, length_cutoff, cpu, dust)
            if best_hit:
                hit_seq = best_hit.split()[-1]
            else:
                print(f'{mirid} not found in query genome')
                continue

            # re-BLAST
            best_rb = run_blast(hit_seq, ref_path, length_cutoff, cpu, dust)
            if best_rb:
                hit_chrom, hit_start, hit_end, hit_strand, hit_length, hit_seq = best_rb.split()
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
            hit.seq.replace('-', '')

            # check location of hit
            print(f'mirna: {mirna.start} - {mirna.end}')
            print(f'hit: {hit.start} - {hit.end}')
            if loc_check(mirna, hit):
                print('CONFIRMED')
                outstr = f'>{mirid}|{hit.chromosome}|{hit.start}|{hit.end}|{hit.strand}\n{hit.seq}\n'
                of.write(outstr)
            else:
                print('REJECTED: Position of best re-BLAST hit does not overlap with reference')


if __name__ == "__main__":
    main()
