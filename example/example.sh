../target/debug/perplexity eval \
    --output bgi-A1.yml \
    --smoothing 1e-8 \
    --test_ecs data/test/eq_classes.txt.gz \
    --train_ecs data/train/eq_classes.txt.gz \
    --quants data/A1-quant.sf