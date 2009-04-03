rc_dict = dict([(f, r) for f, r in zip('ACGT', 'TGCA')])
for f, r in rc_dict.items():
    rc_dict[f.lower()] = r
    
def reverse_complement(dna_string):

    '''Compute the reverse_complement of a string of DNA.  In pure
    python.'''
    
    complement = [rc_dict.get(c, '-') for c in dna_string]
    complement.reverse()
    return ''.join(complement)
