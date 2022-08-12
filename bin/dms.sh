mkdir -p target/log
mkdir -p figures/esm_scatter

# DMS with ESM vote at different stringency cutoffs.
python bin/dms.py > target/log/dms.log
alphas=( 0 0.1 0.2 0.5 0.7 0.9 )
for alpha in ${alphas[@]}
do
    python bin/dms.py --alpha $alpha > target/log/dms_alpha"$alpha".log
done

# Plot stringency experiments.
python bin/plot_dms_alpha_sweep.py

# Compute ESM values for DMS experiments.
python bin/dms_esm.py

# Scatter plots.
python bin/plot_dms_esm_scatter.py

# Compare to Livesey and Marsh across DMS experiments.
python bin/plot_dms_continuous.py
