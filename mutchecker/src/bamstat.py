from yxutil import cmd_run, have_file, rmdir
import numpy as np


def get_mRNA_depth(bam_file, mRNA, depth_file):
    cds_list = [cds for cds in mRNA.sub_features if cds.type == 'CDS']

    if have_file(depth_file):
        rmdir(depth_file)

    for cds in cds_list:
        cmd_string = "samtools depth -Q 30 -aa -r %s:%s-%s %s >> %s" % (
            mRNA.chr_id, cds.start, cds.end, bam_file, depth_file)
        cmd_run(cmd_string)

    return depth_file


def parse_depth_file(depth_file):
    with open(depth_file, 'r') as f:
        depth_list = [int(line.strip().split('\t')[-1]) for line in f]
    depth_list = np.array(depth_list)
    coverage = np.sum(depth_list > 0) / len(depth_list)
    depth = np.mean(depth_list)
    depth_sd = np.std(depth_list)
    return {
        'coverage': coverage,
        'depth': depth,
        'depth_sd': depth_sd
    }


if __name__ == '__main__':
    genome_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.0.fa'
    gff_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Sbicolor.v5.1/Sbicolor_730_v5.1.gene_exons.gff3'
    # bam_file = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/raw_data/map/IFJD.sorted.markdup.bam'
    bam_file = "/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/bam/IPDE.sorted.markdup.bam"
    gene_id = 'Sobic.005G213600.v5.1'
    # gene_id = 'Sobic.001G000200.v5.1'
    work_dir = '/lustre/home/xuyuxing/Work/Jesse/local_adaptation/0.reference/Data/reseq/bam/test'

    from yxseq import read_gff_file, Gene
    from yxutil import cmd_run, have_file, rmdir

    gene_dict = read_gff_file(gff_file)['gene']
    gene = gene_dict[gene_id]
    gene = Gene(from_gf=gene)
    gene.build_gene_seq(genome_file)
    mRNA = gene.model_mRNA

    depth_file = "%s/depth.txt" % work_dir
    get_mRNA_depth(bam_file, mRNA, depth_file)
    parse_depth_file(depth_file)
