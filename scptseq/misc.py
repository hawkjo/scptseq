import os
import re
import gzip
import yaml


ins_re = re.compile('^I(\d+)([ACGT]+)$')
def parse_mutation(mut):
    """
    Return struct of form 
       [mut_type, ref_pos_start, ref_pos_end, bases]
       
    Del has bases == None
    Ins and sub have same ref_pos_start and ref_pos_end
    Sub has tuple of ref and new bases
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


class TargetInfo():
    def __init__(self, target_info_file, gene_name):
        target_info = yaml.load(open(target_info_file), Loader=yaml.FullLoader)
        self.goi_target_info = target_info[gene_name]
        self.gene_chrm, self.gene_start, self.gene_end = [self.goi_target_info[s] for s in ['chrm', 'start', 'end']]
        self.seg_bases = self.goi_target_info['seg_info']['seg_bases']
        self.seg_sites = self.goi_target_info['seg_info']['seg_sites']
        self.skip_sites = self.goi_target_info['seg_info']['skip_sites']
        self.t1_cutsite = self.goi_target_info['targets']['t1']['cutsite']
        self.t2_cutsite = self.goi_target_info['targets']['t2']['cutsite']
        self.cutsites = [self.goi_target_info['targets'][tname]['cutsite'] for tname in ['t1', 't2']]
        self.sorted_seg_sites = sorted(self.seg_sites, key=lambda pos: abs(pos - self.t1_cutsite))
        self.seg_bases_given_site = {
                seg_site: site_seg_bases for seg_site, site_seg_bases in zip(self.seg_sites, self.seg_bases)
                }

    def mut_min_dist_to_cutsite(self, mut):
        mut_type, start, end, bases = parse_mutation(mut)
        if mut_type == 'splice' and isinstance(bases, (int, int)):
            start += bases[0]
            end += bases[1]
        return min([abs(pos - cutsite) for pos in range(start, end+1) for cutsite in self.cutsites])

    def maternal_or_paternal(self, seg_site_bases):
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
