def print_if_full(region):
    if region['type'] is not None and \
       region['gene'] is not None and \
       region['allele'] is not None and \
       region['translation'] is not None:
        fields = [
            region['type'],
            region['gene'],
            region['allele'],
            region['translation'],
        ]
        print('\t'.join([ str(field) for field in fields ]))
        region['type'] = region['gene'] = region['allele'] = region['translation'] = None


def read_field(line, f):
    if line.count('"') == 2:
        return fields[-1].split('"')[-2]
    elif line.count('"') == 0:
        line = f.readline()

    text_field = fields[-1].split('"')[-1].strip()
    for line in f:
        assert(line.startswith('FT'))
        text = line[2:].strip()
        if text.endswith('"'):
            text_field += text.rstrip('"')
            break
        text_field += text

    return text_field
        
if __name__ == '__main__':
    region = {
        'type': None,
        'gene': None,
        'allele': None,
        'translation': None,
    }
    with open('data/imgt/imgt.dat') as f:
        for line in f:
            if not line.startswith('FT'):
                continue
            fields = line.rstrip().split()
            if len(fields) == 3:
                if fields[1] == 'V-REGION':
                    region['type'] = 'V'
                    region['gene'] = region['allele'] = region['translation'] = None
                elif fields[1] == 'D-REGION':
                    region['type'] = 'D'
                    region['gene'] = region['allele'] = region['translation'] = None
                elif fields[1] == 'J-REGION':
                    region['type'] = 'J'
                    region['gene'] = region['allele'] = region['translation'] = None
            elif len(fields) == 2:
                if '/translation' in line:
                    region['translation'] = read_field(line, f)
                    print_if_full(region)
                elif '/gene' in line or '/IMGT_gene' in line:
                    region['gene'] = read_field(line, f)
                    print_if_full(region)
                elif '/allele' in line or '/IMGT_allele' in line:
                    region['allele'] = read_field(line, f)
                    print_if_full(region)
