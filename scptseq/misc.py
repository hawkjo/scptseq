import os
import re
import gzip
import yaml


ins_re = re.compile(r'^I(\d+)([ACGT]+)$')
def parse_mutation(mut: str) -> list:
    """Parse a mutation string into its components.

    Definitional for the mutation notation used throughout the package.
    All positions are 0-based reference coordinates.

    Args:
        mut: A mutation string in one of five forms:

            * {ref}{pos}{alt} — substitution, e.g. G20995489A.
            * D{start}-{end} — deletion, both ends inclusive.
            * I{pos}{bases} — insertion of `bases` immediately after `pos`.
            * J{start}_{end} — standard splice junction, where the bounds are
              the first and last intronic bases.
            * J{start}_{end}>{delta_start}_{delta_end} — a junction offset from
              the standard junction at `(start, end)`. Adding the deltas gives the
              observed intron. J{start}_{end}>- instead means that standard
              junction is absent.

    Returns:
        [mut_type, ref_pos_start, ref_pos_end, bases], where

        `mut_type` is one of 'sub', 'del', 'ins' or 'splice'
        `bases` depends on the type:
            `(ref_base, alt_base)` for a substitution
            `None` for a deletion or standard junction
            The inserted bases for an insertion
            `(delta_start, delta_end)` for an offset junction
            The string 'missing' for an absent junction.
        Substitutions and insertions have `ref_pos_start == ref_pos_end`.

        For an offset junction the positions returned are those of the *standard*
        junction, not of the observed intron.

    Any string not beginning with J, D or I is parsed as a substitution.
    """
    if mut.startswith('J'):
        if '>' not in mut:  # Standard splice location
            us_pos = mut.index('_')
            return ['splice', int(mut[1:us_pos]), int(mut[us_pos+1:]), None]
        arrow_pos = mut.index('>')
        pos_str, delta_str = mut[1:arrow_pos], mut[arrow_pos+1:]
        us_pos = pos_str.index('_')
        start, end = int(pos_str[:us_pos]), int(pos_str[us_pos+1:])
        if delta_str == '-':  # Splice missing
            return ['splice', start, end, 'missing']
        else:  # Splicing error
            us_pos = delta_str.index('_')
            return ['splice', start, end, (int(delta_str[:us_pos]), int(delta_str[us_pos+1:]))]
    elif mut.startswith('D'):
        dash_pos = mut.index('-')
        return ['del', int(mut[1:dash_pos]), int(mut[dash_pos+1:]), None]
    elif mut.startswith('I'):
        m = ins_re.match(mut)
        pos = int(m.group(1))
        return ['ins', pos, pos, m.group(2)]
    else:
        pos = int(mut[1:-1])
        return ['sub', pos, pos, (mut[0], mut[-1])]


def bc_from_fpath(fpath):
    bname = os.path.basename(fpath)
    return bname[:bname.index('.')]


def gzip_friendly_open(fpath, mode='rt'):
    if fpath.endswith('.gz'):
        return gzip.open(fpath, mode)
    return open(fpath, mode)


target_name_re = re.compile(r'^t[1-9]\d*$')
class TargetInfo():
    """One gene's entry from a target-info YAML file.

    These values are read from the file: `chrm`, `start`, `end`, `seg_info.seg_sites`,
    `seg_info.seg_bases`, `seg_info.skip_sites`, and the `cutsite` of each of the targets. Any
    other field present in the file is ignored. All coordinates are 0-based; `start` and `end`
    describe a half-open interval. See `docs/inputs.md` for the full format.

    A gene may define any number of targets, which must be named `t1`, `t2`, ... `tN`. `t1`
    is always required as it is the reference point for ordering the segregating sites.

    Attributes:
        gene_chrm: Contig name, which must match both the BAM header and the genome
            FASTA record id.
        gene_start: Start of the analysed region, inclusive.
        gene_end: End of the analysed region, exclusive.
        seg_sites: Segregating-site positions.
        seg_bases: `[maternal, paternal]` base pair for each segregating site,
            paired with `seg_sites` by list position, in matching order.
        skip_sites: Individual reference positions to ignore when collecting
            mutations.
        target_names: Target names ordered by their numeric suffix, `['t1', 't2', ...]`
        cutsites: The cut site of each target, in `target_names` order. 
        t1_cutsite: Cut site of `t1`, the reference target.
        sorted_seg_sites: Segregating sites ordered by distance from the `t1` cut
            site, which is the order `maternal_or_paternal` consults them in.
        seg_bases_given_site: Maps each segregating site to its base pair.

    Raises:
        ValueError: If the gene's targets are not named `t1`, `t2`, ... `tN`, or if `t1`
            is absent.
    """
    def __init__(self, target_info_file, gene_name):
        target_info = yaml.safe_load(open(target_info_file))
        self.goi_target_info = target_info[gene_name]
        self.gene_chrm, self.gene_start, self.gene_end = [self.goi_target_info[s] for s in ['chrm', 'start', 'end']]
        self.seg_bases = self.goi_target_info['seg_info']['seg_bases']
        self.seg_sites = self.goi_target_info['seg_info']['seg_sites']
        self.skip_sites = self.goi_target_info['seg_info']['skip_sites']
        targets = self.goi_target_info['targets']
        if 't1' not in targets or not all(target_name_re.match(str(tname)) for tname in targets):
            raise ValueError(
                    f'{gene_name}: targets must be named t1, t2, ... with t1 present. '
                    f'Got: {sorted(map(str, targets))}'
                    )
        # Ordered by target number rather than lexicographic or by order in the file
        self.target_names = sorted(targets, key=lambda tname: int(tname[1:]))
        self.cutsites = [targets[tname]['cutsite'] for tname in self.target_names]
        self.t1_cutsite = targets['t1']['cutsite']
        self.sorted_seg_sites = sorted(self.seg_sites, key=lambda pos: abs(pos - self.t1_cutsite))
        self.seg_bases_given_site = {
                seg_site: site_seg_bases for seg_site, site_seg_bases in zip(self.seg_sites, self.seg_bases)
                }

    def mut_min_dist_to_cutsite(self, mut: str) -> int:
        """Distance in bases from a mutation to the nearest of the cut sites.

        Args:
            mut: A mutation string, as parsed by `parse_mutation`.

        Returns:
            The smallest distance in bases to any cut site.
        """
        mut_type, start, end, bases = parse_mutation(mut)
        if mut_type == 'splice' and isinstance(bases, (int, int)):
            start += bases[0]
            end += bases[1]
        return min([abs(pos - cutsite) for pos in range(start, end+1) for cutsite in self.cutsites])

    def maternal_or_paternal(self, seg_site_bases: dict):
        """Assign a haplotype from the bases a read shows at the segregating sites.

        Segregating sites are consulted in order of distance from the `t1` cut site. The first site
        the read covers whose base matches one of the two known alleles gives the haplotype. A site
        covered but matching neither allele falls through to the next site.

        Args:
            seg_site_bases: Maps segregating-site position to the base observed at
                that position in one read. Sites the read does not cover are absent.

        Returns:
            'maternal', 'paternal', or `None` if no covered site matches either allele

        Which allele is maternal is supplied in the target file, not derived from the reference
        sequence.
        """
        for seg_site in self.sorted_seg_sites:
            if seg_site not in seg_site_bases:
                continue
            obs_base = seg_site_bases[seg_site]
            site_bases = self.seg_bases_given_site[seg_site]
            if obs_base == site_bases[0]:
                return 'maternal'
            if obs_base == site_bases[1]:
                return 'paternal'
        return None
