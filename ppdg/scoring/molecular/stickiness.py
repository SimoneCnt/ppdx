
def get_stickiness(res, klass):
    """
        Compute the average stickiness of the core interface residues according
        to Levy [1].
        [1] E. D. Levy, S. De, and S. A. Teichmann, "Cellular crowding imposes 
            global constraints on the chemistry and evolution of proteomes", 
            PNAS, vol. 109, no. 50, pp. 20461-20466, 2012.
    """
    stk = {
        'ALA': 0.0062, 'CYS': 1.0372, 'ASP':-0.7485, 'GLU':-0.7893, 'PHE': 1.2727,
        'GLY':-0.1771, 'HIS': 0.1204, 'ILE': 1.1109, 'LYS':-1.1806, 'LEU': 0.9138,
        'MET': 1.0124, 'ASN':-0.2693, 'PRO':-0.1799, 'GLN':-0.4114, 'ARG':-0.0876,
        'SER': 0.1376, 'THR': 0.1031, 'VAL': 0.7599, 'TRP': 0.7925, 'TYR': 0.8806
        }
    ll = list()
    for r, k in zip(res, klass):
        if k=='core':
            ll.append(stk[r])
    return np.mean(ll)

