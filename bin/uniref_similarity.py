from Bio import SeqIO
import sys
from thefuzz import fuzz

if __name__ == '__main__':
    fname_database = sys.argv[1]
    seq_target = sys.argv[2]

    closest = {}
    min_score = -1
    n_closest, max_closest = 0, 10000
    
    for record in SeqIO.parse(fname_database, 'fasta'):
        seq = str(record.seq)
        similarity_score = fuzz.ratio(seq_target, seq)
        
        # If collection is too small or
        # if sequence is better than collection members,
        # then add sequence to collection.
        if n_closest < max_closest or \
           similarity_score >= min_score:
            if similarity_score not in closest:
                closest[similarity_score] = []
            closest[similarity_score].append(seq)
            n_closest += 1

            if min_score == -1: # Edge case at beginning.
                min_score = similarity_score
            elif n_closest <= max_closest and similarity_score < min_score:
                min_score = similarity_score
                                          
        # If collection gets too big, pop lowest ones off.
        while n_closest > max_closest and len(closest) > 1:
            n_closest -= len(closest.pop(min_score))
            assert(n_closest >= 0)
            min_score = min(closest.keys())

    # Output closest sequences.
    for similarity_score in reversed(sorted(closest.keys())):
        fields = [ similarity_score, closest[similarity_score] ]
        print('\t'.join([ str(field) for field in fields ]))
