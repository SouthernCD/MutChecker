from mutchecker.src.utils import read_sequence_from_file
from mutchecker.src.assem import get_range_haplotype, get_all_hap_combinations
from mutchecker.src.miniprot import miniprot_job
from mutchecker.src.bamstat import get_mRNA_depth, parse_depth_file, get_all_cds_depth
from yxseq import read_gff_file, Gene
from yxmath.interval import merge_intervals, overturn, interval_minus_set
from yxutil import cmd_run, mkdir, rmdir, pickle_load, have_file, multiprocess_running
import json
import os
import shutil
import warnings
warnings.filterwarnings('ignore')


def reseq_pipeline(gene_id, mRNA, genome_file, gff_file, bam_file, exon_extend, work_dir, clean):

    # skip if the gene has been processed
    if not have_file(work_dir + "/result.json"):
        # Load gene information
        mkdir(work_dir, False)
        cds_list = [cds for cds in mRNA.sub_features if cds.type == 'CDS']

        # state reseq depth and coverage in the CDS region
        depth_file = "%s/cds_depth.txt" % work_dir
        get_mRNA_depth(bam_file, mRNA, depth_file)
        depth_report = parse_depth_file(depth_file)

        # assemble gene sequences
        ass_work_dir = work_dir + "/assembly"
        mkdir(ass_work_dir)

        # assemble each range
        range_list = interval_minus_set((mRNA.start, mRNA.end), overturn(merge_intervals(
            [(cds.start-exon_extend, cds.end+exon_extend) for cds in cds_list])))
        # range_list = merge_intervals(
        #     [(cds.start-exon_extend, cds.end+exon_extend) for cds in cds_list])

        hap_seq_dict = {}
        ref_seq_dict = {}
        for i, (start, end) in enumerate(range_list):
            range_dir = "%s/range%d" % (ass_work_dir, i+1)
            mkdir(range_dir)
            chr_id = mRNA.chr_id
            hap1_file, hap2_file, ref_file = get_range_haplotype(
                i, chr_id, start, end, bam_file, genome_file, range_dir)
            hap_seq_dict[i] = (hap1_file, hap2_file)
            ref_seq_dict[i] = ref_file

        all_combination_sequences = list(
            set(get_all_hap_combinations(hap_seq_dict)))
        ref_sequence = ''.join([read_sequence_from_file(ref_file)
                                for ref_file in ref_seq_dict.values()])

        hap_gene_dict = {}
        for i, seq in enumerate(all_combination_sequences):
            hap_gene_dict[i] = "%s/haplotype%d_genome.fasta" % (
                ass_work_dir, i+1)
            with open(hap_gene_dict[i], 'w') as f:
                f.write(f">haplotype{i+1}\n{seq}\n")

        ref_seq_file = "%s/reference_genome.fasta" % ass_work_dir
        with open(ref_seq_file, 'w') as f:
            f.write(f">ref_genome\n{ref_sequence}\n")

        # clustalw nucleotide alignment
        genome_aln_file = "%s/all_genome.fasta" % ass_work_dir
        with open(genome_aln_file, 'w') as f:
            f.write(f">ref_genome\n{ref_sequence}\n")
            for i, seq in enumerate(all_combination_sequences):
                f.write(f">haplotype{i+1}\n{seq}\n")

        cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=DNA" % (
            genome_aln_file, genome_aln_file)
        cmd_run(cmd_string, cwd=ass_work_dir)

        # re-annotation
        anno_work_dir = work_dir + "/annotation"
        mkdir(anno_work_dir)

        # write reference query protein sequence
        ref_prot_file = "%s/reference_protein.fasta" % anno_work_dir
        with open(ref_prot_file, 'w') as f:
            f.write(">%s\n%s\n" % (gene_id, mRNA.aa_seq))

        ref_cds_file = "%s/reference_CDS.fasta" % anno_work_dir
        with open(ref_cds_file, 'w') as f:
            f.write(">%s\n%s\n" % (gene_id, mRNA.cds_seq))

        # run miniprot
        reanno_ref_pseudo_dict, reanno_ref_pt_seq, reanno_ref_cds_seq = miniprot_job(
            ref_prot_file, ref_seq_file, anno_work_dir + "/reference_miniprot")

        reanno_ref_pt_file = "%s/reanno_ref_protein.fasta" % anno_work_dir
        with open(reanno_ref_pt_file, 'w') as f:
            f.write(">reanno_ref\n%s\n" % reanno_ref_pt_seq)

        reanno_ref_cds_file = "%s/reanno_ref_CDS.fasta" % anno_work_dir
        with open(reanno_ref_cds_file, 'w') as f:
            f.write(">reanno_ref\n%s\n" % reanno_ref_cds_seq)

        hap_miniport_results_dict = {}
        for i, hap_gene_file in hap_gene_dict.items():
            pseudo_dict, pt_seq, cds_seq = miniprot_job(
                ref_prot_file, hap_gene_file, anno_work_dir + f"/haplotype{i+1}_miniprot")
            hap_miniport_results_dict[i] = (pseudo_dict, pt_seq, cds_seq)

            hap_pt_file = anno_work_dir + f"/haplotype{i+1}_protein.fasta"
            with open(hap_pt_file, 'w') as f:
                f.write(f">haplotype{i+1}\n{pt_seq}\n")

            hap_cds_file = anno_work_dir + f"/haplotype{i+1}_CDS.fasta"
            with open(hap_cds_file, 'w') as f:
                f.write(f">haplotype{i+1}\n{cds_seq}\n")

        # run clustalw protein alignment
        all_protein_file = anno_work_dir + "/all_protein.fasta"
        with open(all_protein_file, 'w') as f:
            f.write(">%s\n%s\n" % (gene_id, mRNA.aa_seq))
            f.write(f">reanno_ref\n{reanno_ref_pt_seq}\n")
            for i, (_, pt_seq, _) in hap_miniport_results_dict.items():
                f.write(f">haplotype{i+1}\n{pt_seq}\n")

        # cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=PROTEIN" % (
        #     all_protein_file, all_protein_file)
        cmd_string = "mafft --preservecase --auto %s > %s.aln" % (
            all_protein_file, all_protein_file)
        cmd_run(cmd_string, cwd=anno_work_dir)

        # run clustalw cds alignment
        all_cds_file = anno_work_dir + "/all_CDS.fasta"
        with open(all_cds_file, 'w') as f:
            f.write(">%s\n%s\n" % (gene_id, mRNA.cds_seq))
            f.write(
                f">reanno_ref\n{read_sequence_from_file(reanno_ref_cds_file)}\n")
            for i, (_, _, cds_seq) in hap_miniport_results_dict.items():
                f.write(f">haplotype{i+1}\n{cds_seq}\n")

        # cmd_string = "clustalw2 -INFILE=%s -ALIGN -OUTPUT=FASTA -OUTFILE=%s.aln -type=DNA" % (
        #     all_cds_file, all_cds_file)
        cmd_string = "mafft --preservecase --auto %s > %s.aln" % (
            all_cds_file, all_cds_file)
        cmd_run(cmd_string, cwd=anno_work_dir)

        # merge all results
        results_dict = {'gene_id': gene_id, 'genome_file': genome_file, 'gff_file': gff_file,
                        'bam_file': bam_file, 'exon_extend': exon_extend}
        results_dict['reseq_stat'] = depth_report
        results_dict['ref_redo'] = reanno_ref_pseudo_dict

        for i in hap_miniport_results_dict:
            results_dict[f"haplotype{i+1}"] = hap_miniport_results_dict[i][0]

        low_coverage = True if results_dict['reseq_stat']['coverage'] < 0.5 else False
        low_depth = True if results_dict['reseq_stat']['depth'] < 5 else False

        fatal_mut = True
        for i in results_dict:
            if i.startswith('haplotype'):
                if results_dict[i]['frameshift'] is False and results_dict[i]['stopcodon'] is False and results_dict[i]['headmissed'] is False and results_dict[i]['coverage'] > 0.9 and results_dict[i]['identity'] > 0.9:
                    fatal_mut = False
                    break

        ref_bad_mut = results_dict['ref_redo']['frameshift'] or results_dict['ref_redo'][
            'stopcodon'] or results_dict['ref_redo']['headmissed'] or results_dict['ref_redo']['coverage'] < 0.9
        fatal_mut = False if ref_bad_mut else fatal_mut
        fatal_mut = False if low_depth else fatal_mut

        if low_coverage or fatal_mut:
            results_dict['dead'] = True
        else:
            results_dict['dead'] = False

        result_json = "%s/result.json" % work_dir
        with open(result_json, 'w') as f:
            f.write(json.dumps(results_dict, indent=4))

    else:
        print("Gene %s has been processed, skip." % gene_id)

    if clean:
        # assembly
        for i in os.listdir("%s/assembly" % work_dir):
            if not i.startswith('all'):
                rmdir("%s/assembly/%s" % (work_dir, i))
            elif i.endswith('.dnd'):
                rmdir("%s/assembly/%s" % (work_dir, i))
        # annotation
        for i in os.listdir("%s/annotation" % work_dir):
            if not i.startswith('all') and not i.endswith('_miniprot'):
                rmdir("%s/annotation/%s" % (work_dir, i))
            elif i.endswith('.dnd'):
                rmdir("%s/annotation/%s" % (work_dir, i))
            elif i.endswith('_miniprot'):
                rmdir("%s/annotation/%s/miniprot.out" % (work_dir, i))
                rmdir("%s/annotation/%s/anno_pt.fasta" % (work_dir, i))
        # depth
        rmdir("%s/cds_depth.txt" % work_dir)

        # zip work_dir
        shutil.make_archive(work_dir, 'zip', work_dir)
        rmdir(work_dir)


def reseq_main(args):
    genome_file = os.path.realpath(args.genome_file)
    gff_file = os.path.realpath(args.gff_file)
    bam_file = os.path.realpath(args.bam_file)
    gene_id = args.gene_id
    work_dir = os.path.realpath(args.output_dir)
    exon_extend = args.exon_extend
    clean = args.clean

    if not gene_id is None:
        if os.path.isfile(gff_file):
            gene_dict = read_gff_file(gff_file)['gene']
            gene = gene_dict[gene_id]
            gene = Gene(from_gf=gene)
            gene.build_gene_seq(genome_file)
            mRNA = gene.model_mRNA
        elif os.path.isdir(gff_file):
            mRNA_pkl = "%s/%s.pkl" % (gff_file, gene_id)
            mRNA = pickle_load(mRNA_pkl)
        reseq_pipeline(gene_id, mRNA, genome_file, gff_file,
                       bam_file, exon_extend, work_dir, clean)
    else:
        gene_dict = read_gff_file(gff_file)['gene']
        mkdir(work_dir)
        num = 0
        for gene_id in gene_dict:
            gene = gene_dict[gene_id]
            gene = Gene(from_gf=gene)
            gene.build_gene_seq(genome_file)
            mRNA = gene.model_mRNA
            reseq_pipeline(gene_id, mRNA, genome_file, gff_file,
                           bam_file, exon_extend, work_dir+"/"+gene_id, clean)
            num += 1
            print("Gene %s done (%d/%d)." % (gene_id, num, len(gene_dict)))


def bamstat_main(args):
    depth_file = args.bam_file + ".cds.depth"
    result_json = args.bam_file + ".cds.depth.json"
    if not have_file(depth_file) and not have_file(result_json):
        get_all_cds_depth(args.bam_file, args.gff_file, depth_file)
        depth_output_dict = parse_depth_file(depth_file)
        with open(result_json, 'w') as f:
            f.write(json.dumps(depth_output_dict, indent=4))


if __name__ == '__main__':

    class abc():
        pass

    args = abc()
    args.genome_file = 'Sbicolor_730_v5.0.fa'
    args.gff_file = 'Sbicolor_730_v5.1.gene_exons.gff3'
    args.bam_file = "IZPX.sorted.markdup.bam"
    args.gene_id = 'Sobic.003G269600.v5.1'
    args.output_dir = 'Sobic.003G269600.v5.1.vs.IZPX'
    args.exon_extend = 500
    args.clean = True

    reseq_main(args)
