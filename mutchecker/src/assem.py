from yxutil import cmd_run
import itertools
from mutchecker.src.utils import read_sequence_from_file


def get_range_haplotype(range_id, chr_id, start, end, raw_bam_file, genome_file, work_dir):
    """
    Get the haplotype sequences of a specific region.
    Parameters:
    - range_id: An integer to identify the region
    - chr_id: Chromosome name
    - start: Start position of the region
    - end: End position of the region
    - raw_bam_file: Path to the original BAM file
    - genome_file: Path to the reference genome file
    - work_dir: Path to the working directory
    """

    cmd_string = "samtools view -bS %s %s:%d-%d > range%d.bam" % (
        raw_bam_file, chr_id, start, end, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "samtools index range%d.bam" % (range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "freebayes -f %s range%d.bam > range%d_variants.vcf" % (
        genome_file, range_id, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "whatshap phase -o range%d_phased.vcf --reference=%s range%d_variants.vcf range%d.bam" % (
        range_id, genome_file, range_id, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bgzip range%d_phased.vcf && tabix range%d_phased.vcf.gz" % (
        range_id, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "samtools faidx %s %s:%d-%d > range%d.ref.fa" % (
        genome_file, chr_id, start, end, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bcftools consensus -H 1 -f range%d.ref.fa range%d_phased.vcf.gz > range%d_hap1.fasta" % (
        range_id, range_id, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bcftools consensus -H 2 -f range%d.ref.fa range%d_phased.vcf.gz > range%d_hap2.fasta" % (
        range_id, range_id, range_id)
    cmd_run(cmd_string, cwd=work_dir)

    hap1_file = "%s/range%d_hap1.fasta" % (work_dir, range_id)
    hap2_file = "%s/range%d_hap2.fasta" % (work_dir, range_id)
    ref_file = "%s/range%d.ref.fa" % (work_dir, range_id)

    return hap1_file, hap2_file, ref_file


def get_all_hap_combinations(hap_seq_dict):
    # 获取所有haplotype的组合
    all_combinations = list(itertools.product(*hap_seq_dict.values()))
    # 存储每个组合的序列
    combination_sequences = []

    for combination in all_combinations:
        # 对于每个组合，连接其序列
        combined_sequence = ''.join(
            [read_sequence_from_file(hap_file) for hap_file in combination])
        combination_sequences.append(combined_sequence)

    return combination_sequences
