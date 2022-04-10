declare -a models=( "esm1b" "esm1v1" "esm1v2" "esm1v3" "esm1v4" "esm1v5" )

for model in "${models[@]}"
do
    echo $model
    python bin/$1.py --model-name $model > reconstruct_$model.log
done

bash bin/grep_esm1v.sh .
