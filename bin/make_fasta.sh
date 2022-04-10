line_num=0

cat $1 | \
    while read line
    do
        echo ">seq"$line_num
        echo $line
        ((line_num+=1))
    done
