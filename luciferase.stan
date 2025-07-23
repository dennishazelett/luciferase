data {
    int nrow;
    int nenh;
    int nplas;
    int ntrsf;
    int nbatch;
    int ff[nrow];
    int rn[nrow];
    int enhancer[nrow];
    int transfection[nrow];
    int plasmid[nrow];
    int batch[nrow];
    int dht[nrow];
}
transformed data {
    real ld[nrow];
    int dhtidx[nrow];
    int nenhTimes2;
}
parameters {
    // estimate of the effect of the enhancer
    real enhEff[nenh];
    // estimate dht effect for each enhancer
    real dhtEff[nenh];
    // estimate plasmid preparation effect
    real plpEff[nplas];
    // estimate transfection effect
    real traneff[ntrsf];
    // estimate batch effect
    real batchEff[nbatch];
    // estimate hyperparams on transformation and plasmid prep
    real<lower=0> sdplas;
    real<lower=0> trsd;
    // additional error
    real<lower=0> sd; // on log(firefly/renilla)
    real<lower=0> sdbatch; // st dev of batch
    real<lower=0> sdenh;   // st dev of enhancer
    real<lower=0> sddht;   // st dev of dht effect
    real<lower=2> enhdof;  // enhancer deg freedom
    real trsfdof;          // transfection deg freedom
}
model {
    // prior:
    sdenh ~ exponential(1.0/8.0);
    enhdof ~ exponential(1.0/20.0);
    for (i in 1:nenh) {
           enhEff[i] ~ student_t(enhdof, 0, sdenh);
    }
    // prior on dht effect and its st dev
    sddht ~ exponential(1/0.1);
    for (i in 1:nenh) {
        dhtEff[i] ~ cauchy(0, sddht);
    }
    sdbatch ~ exponential(1.0/1.0);
    for (i in 1:nbatch) {
        batchEff[i] ~ normal(0, sdbatch);
    }
    trsd ~ exponential(1.0/0.5);
    trsfdof ~ gamma(1.1, 3);
    for(i in 1:ntrsf) {
        traneff[i] ~ student_t(trsfdof, 0, trsd);
    }
    sdplas ~ exponential(1.0/1.0);
    for (i in 1:nplas) {
        plpEff[i] ~ normal(0, sdplas);
    }
    sd ~ exponential(1.0/1.0);
    // likelihood fxn:
    for (i in 1:nrow) {
        ld[i] ~ normal(mean(ld) + enhEff[enhancer[i]] + (dht[i] - 1) * dhtEff[enhancer[i]] + plpEff[plasmid[i]] + traneff[transfection[i]] + batchEff[batch[i]], sd);
    }
}
