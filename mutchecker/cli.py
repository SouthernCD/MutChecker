import argparse
from mutchecker.src.pipeline import reseq_main, bamstat_main


class CustomHelpFormatter(argparse.HelpFormatter):
    def add_subparsers(self, *args, **kwargs):
        subparsers_action = super().add_subparsers(*args, **kwargs)
        subparsers_action._parser_class = CustomSubcommandParser
        return subparsers_action


class CustomSubcommandParser(argparse.ArgumentParser):
    def format_help(self):
        formatter = self._get_formatter()

        # Add the usage
        formatter.add_usage(self.usage, self._actions,
                            self._mutually_exclusive_groups)

        # Add the description
        formatter.add_text(self.description)

        # Add the subcommands
        for action in self._actions:
            if isinstance(action, argparse._SubParsersAction):
                formatter.start_section("subcommands")
                for choice, subparser in action.choices.items():
                    formatter.add_text(f"{choice}: {subparser.description}\n")
                formatter.end_section()

        # Add the epilog
        formatter.add_text(self.epilog)

        # Return the full help string
        return formatter.format_help()


class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse

        parser = argparse.ArgumentParser(
            prog='mutchecker',
            description="Main command description.",
            formatter_class=CustomHelpFormatter
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for bamstat
        parser_a = subparsers.add_parser('bamstat',
                                         description='Stat the depth and coverage of all CDS in a bam file',
                                         help='Stat the depth and coverage of all CDS in a bam file')
        parser_a.add_argument('gff_file', type=str,
                              help='reference gff file')
        parser_a.add_argument('bam_file', type=str,
                              help='sorted and markdup bam file')

        parser_a.set_defaults(func=bamstat_main)

        # argparse for reseq
        parser_a = subparsers.add_parser('reseq',
                                         description='Check the mutation of a specific gene by resequencing data',
                                         help='Check the mutation of a specific gene by resequencing data')
        parser_a.add_argument('genome_file', type=str,
                              help='reference genome file')
        parser_a.add_argument('gff_file', type=str,
                              help='reference gff file')
        parser_a.add_argument('bam_file', type=str,
                              help='sorted and markdup bam file')
        parser_a.add_argument('-i', '--gene_id', type=str,
                              help='gene id, default=None', default=None)
        parser_a.add_argument('-o', '--output_dir', type=str,
                              help='output directory, default=\"mutchecker_output\"',
                              default="mutchecker_output")
        parser_a.add_argument('-e', '--exon_extend', type=int,
                              help='extend length of cds, default=500',
                              default=500)
        parser_a.add_argument('-c', '--clean', action='store_true',
                              help='clean the intermediate files, default False')

        parser_a.set_defaults(func=reseq_main)

        self.arg_parser = parser

        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()

        if self.args.subcommand_name == 'reseq':
            reseq_main(self.args)
        elif self.args.subcommand_name == 'bamstat':
            bamstat_main(self.args)
        else:
            self.arg_parser.print_help()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
