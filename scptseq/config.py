import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return []

    def _expanduser_or_none(self, arg):
        return os.path.expanduser(self._arguments[arg]) if self._arguments[arg] else None

    @property
    def command(self):
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('count', 'splitfastqs'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def fastq_files(self):
        return [os.path.expanduser(fp) for fp in self._comma_delimited_arg('<fastq_files>')]

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def gene_name(self):
        return self._arguments['<gene_name>']

    @property
    def run_name(self):
        return self._arguments['<run_name>']

    @property
    def genome_file(self):
        return self._expanduser_or_none('<genome_file>')

    @property
    def control_bam_dir(self):
        return self._expanduser_or_none('<control_bam_dir>')

    @property
    def perturbed_bam_dir(self):
        return self._expanduser_or_none('<perturbed_bam_dir>')

    @property
    def target_info_file(self):
        return self._expanduser_or_none('<target_info_file>')

    @property
    def threads(self):
        return int(self._arguments['--threads'] or 1)

    @property
    def output_dir(self):
        return self._arguments['--output-dir'] or '.'
