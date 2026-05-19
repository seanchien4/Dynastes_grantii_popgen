#!/bin/bash
#SBATCH --job-name=M0
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --output=M0.out
#SBATCH --error=M0.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schien@tamu.edu
#SBATCH --account=132704493814

# Variables

MODEL="M0_isolation"
SFS_PREFIX="D_grantii"
NREPS=50
NITER=40
NSIMS=10

echo "=== Running Model 0: Isolation ==="
echo "Started: $(date)"

for REP in $(seq 1 $NREPS)
do
    echo "--- Running replicate $REP of $NREPS ---"

    mkdir -p rep_${REP}
    cd rep_${REP}

    cp ../${MODEL}.tpl .
    cp ../${MODEL}.est .
    cp ../${SFS_PREFIX}_MSFS.obs .

    $FSC \
        --tplfile ${MODEL}.tpl \
        --estfile ${MODEL}.est \
        --multiSFS \
        -E 1 \
        --numsims $NSIMS \
        --maxiter $NITER \
        --minSFSCount 1 \
    --cores $SLURM_CPUS_PER_TASK \
        --quiet

    if [ $? -eq 0 ]; then
        echo "Rep $REP: SUCCESS"
        LIKE=$(awk 'NR==2{print $NF}' \
            ${MODEL}/${MODEL}.bestlhoods 2>/dev/null)
        echo "Rep $REP likelihood: $LIKE"
    else
        echo "Rep $REP: FAILED"
    fi

    cd ..
done

echo ""
echo "=== Finding best run ==="
BEST_LIKE=-999999999
BEST_REP=0

for REP in $(seq 1 $NREPS)
do
    LHOOD_FILE="rep_${REP}/${MODEL}/${MODEL}.bestlhoods"
    if [ -f "$LHOOD_FILE" ]; then
        LIKE=$(awk 'NR==2{print $NF}' $LHOOD_FILE)
        IS_BETTER=$(echo "$LIKE > $BEST_LIKE" | bc -l 2>/dev/null)
        if [ "$IS_BETTER" = "1" ]; then
            BEST_LIKE=$LIKE
            BEST_REP=$REP
        fi
    fi
done

echo "Best replicate: $BEST_REP"
echo "Best likelihood: $BEST_LIKE"

if [ $BEST_REP -gt 0 ]; then
    mkdir -p best_run
    cp -r rep_${BEST_REP}/${MODEL}/* best_run/
    echo ""
    echo "=== Best run results ==="
    cat best_run/${MODEL}.bestlhoods
fi

echo "Ended: $(date)"