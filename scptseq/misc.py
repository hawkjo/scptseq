import re


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


