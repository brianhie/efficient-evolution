python bin/dms.py > target/log/dms.log

alphas=( 0 0.1 0.2 0.5 0.7 0.9 )

for alpha in ${alphas[@]}:
do
    python bin/dms.py --alpha $alpha > target/log/dms_alpha"$alpha".log
done

python bin/plot_dms_alpha_sweep.py
