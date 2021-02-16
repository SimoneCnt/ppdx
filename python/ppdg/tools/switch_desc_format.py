
def switch_desc_format(descriptors):
    """
        Convert dict from 
            index - descriptor - value
        to
            descriptor - index - value
        and vice versa.
        For internal use only...
    """
    desc_flat = dict()
    for key, values in descriptors.items():
        for index, value in values.items():
            if not index in desc_flat.keys():
                desc_flat[index] = dict()
            desc_flat[index][key] = float(value)
    return desc_flat

