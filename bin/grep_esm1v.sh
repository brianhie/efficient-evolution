grep -h -P "VH mut" $1/reconstruct_* | \
    cut -f2- -d':' | \
    sed 's/,/\n/g' | \
    sed 's/ //g' | \
    sort | \
    uniq -c | \
    sort -k1,1nr
grep -h -P "VH mut" $1/reconstruct_* | \
    cut -f2- -d':'

grep -h -P "VL mut|VK mut" $1/reconstruct_* | \
    cut -f2- -d':' | \
    sed 's/,/\n/g' | \
    sed 's/ //g' | \
    sort | \
    uniq -c | \
    sort -k1,1nr
grep -h -P "VL mut|VK mut" $1/reconstruct_* | \
    cut -f2- -d':' 
