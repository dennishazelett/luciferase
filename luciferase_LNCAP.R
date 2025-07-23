library(ggplot2)
library(rstan)
library(parallel)
library(grid)
source('~/projects/ProgramsDoingBayesianDataAnalysis/HDIofMCMC.R')

ffdata <- read.csv(file="luciferase_assays.csv", sep=",", header=TRUE)
a <- ffdata
ld <- log(a$firefly/a$renilla)
controlmean <- with(subset(subset(a, enhancer==c("ct1", "ct2")), DHT==1),
                    mean(log(firefly/renilla)))


fflist = list(nrow=dim(a)[1],
    nenh = length(unique(a$enhancer)),
    nplas = length(unique(a$plasmid)),
    ntrsf = length(unique(a$transfection)),
    nbatch = length(unique(a$batch)),
    controlmean = controlmean,
    ff=a$firefly,rn=a$renilla,
    enhancer=as.numeric(a$enhancer),
    transfection=a$transfection,
    plasmid=a$plasmid,
    dht=a$DHT,
    batch=a$batch)



firefly_model <- "
data {
    int nrow;
    int nenh;
    int nplas;
    int ntrsf;
    int nbatch;
    real controlmean;
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
    for (i in 1:nrow) {
        ld[i] <- log(ff[i]) - log(rn[i]);
    }
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
    real<lower=0> trsfdof;          // transfection deg freedom
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
    // likelihood:
    for (i in 1:nrow) {
        ld[i] ~ normal(controlmean + enhEff[enhancer[i]] + (dht[i] - 1) * dhtEff[enhancer[i]] + plpEff[plasmid[i]] + traneff[transfection[i]] + batchEff[batch[i]], sd);
    }
}
"

seed <- 11041974

f1 <- stan(model_code = firefly_model, data = fflist, iter=20000, warmup=20000, thin=1, chains=1, seed=seed, chain_id = 1)

sflist1 <- mclapply(1:4, mc.cores = 8,
                   function(i) stan(fit = f1,
                                    data = fflist,
                                    seed = seed,
                                    iter=200000,
                                    thin=40,
                                    chains = 1,
                                    chain_id = i,
                                    refresh = -1))

fit.a <- sflist2stanfit(sflist1)



foo <- extract(fit.a)
save.image(fit.a, file="fireflysam.R")

pdf("control1plusminusDHT.pdf")
hdi <- HDIofMCMC(foo$dhtEff[,1])
ggplot(data.frame(enhdiff=foo$dhtEff[,1])) + geom_histogram(aes(enhdiff)) + ggtitle("Control 1, DHT effect") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
dev.off()


pdf("control2plusminusDHT.pdf")
hdi <- HDIofMCMC(foo$dhtEff[,2])
ggplot(data.frame(enhdiff=foo$dhtEff[,2])) + geom_histogram(aes(enhdiff)) + ggtitle("Control 2, DHT effect") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
dev.off()


pdf("F35snp.pdf")
hdi <- HDIofMCMC(foo$enhEff[,26] - foo$enhEff[,25])
ggplot(data.frame(enhdiff=foo$enhEff[,26] - foo$enhEff[,25])) + geom_histogram(aes(enhdiff)) + labs(title="F35m - F35 allele comparison, unstimulated") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x = hdi[2], xend = hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))# m35-3 - 35-3
dev.off()

pdf("F30-1plusminusDHT.pdf")
hdi <- HDIofMCMC(foo$dhtEff[,6])
ggplot(data.frame(enhdiff=foo$dhtEff[,6])) + geom_histogram(aes(enhdiff)) + labs(title="F30-1 (rs10455782) IGF2R locus, DHT vs. unstimulated", y="log fold change") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
dev.off()

pdf("F37snp.pdf")
hdi <- HDIofMCMC(foo$enhEff[,32] - foo$enhEff[,31])
ggplot(data.frame(enhdiff=foo$enhEff[,32] - foo$enhEff[,31])) + geom_histogram(aes(enhdiff))  + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3))) # m37-1 - 37-1
dev.off()





################################
# SECTION: ANALYSIS of F38 SNP #
################################

pdf("F38snp.pdf")
hdi <- HDIofMCMC(foo$enhEff[,38] - foo$enhEff[,37])
ggplot(data.frame(enhdiff=foo$enhEff[,38] - foo$enhEff[,37])) + geom_histogram(aes(enhdiff)) + ggtitle("F38m - F38 allele contrast (unstimulated)") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3))) # m38 - 38
dev.off()

hdi <- HDIofMCMC(foo$dhtEff[,38] - foo$dhtEff[,37], cred=0.95)
f38 <- ggplot(data.frame(enhdiff=foo$dhtEff[,38] - foo$dhtEff[,37])) + geom_histogram(aes(enhdiff)) + ggtitle("F38m - F38 allele contrast (DHT effect)") 
f38 <- f38 + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3)
f38 <- f38 + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
pdf("F38snpPlusDHT.pdf")
f38
dev.off()

pdf("F38mPlusMinusDHT.pdf")
stim <- foo$dhtEff[,38]
hdi <- HDIofMCMC(foo$dhtEff[,38])
ggplot(data.frame(enhdiff=stim)) + geom_histogram(aes(enhdiff)) + labs(title="F38m DHT effect") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3))) # m38 +-DHT
dev.off()

stim <- foo$dhtEff[,37]
hdi <- HDIofMCMC(foo$dhtEff[,37])
pdf("F38PlusMinusDHT.pdf") 
ggplot(data.frame(enhdiff=stim)) + geom_histogram(aes(enhdiff)) + labs(title="F38 DHT effect", x="fold change", y="probability mass") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3))) # 38 +-DHT
dev.off()

stimdiff <- foo$dhtEff[,38]-foo$dhtEff[,37]
hdi <- HDIofMCMC(stimdiff)
pdf("F38allelediffDHT.pdf") 
ggplot(data.frame(enhdiff=stimdiff)) + geom_histogram(aes(enhdiff)) + labs(title="F38, effect of rs4907792 on DHT response", x="log fold change", y="probability mass") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3))) # 38 +-DHT
dev.off()



################################
# SECTION: ANALYSIS of PSA     #
################################

stim <- foo$dhtEff[,51]
hdi <- HDIofMCMC(stim)
mstim <- mean(stim)
psa <- ggplot(data.frame(enhdiff=stim))
psa <- psa + geom_histogram(aes(enhdiff)) 
psa <- psa + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=2 )
psa <- psa + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=2 ) 
psa <- psa + geom_text(aes(x=hdi, y=rep(550,2), label=formatC(hdi, digits=3)))
psa <- psa + geom_segment(aes(x=hdi[1], xend=hdi[2], y=100, yend=100), arrow=arrow(length=unit(0.1, "in"), ends="both"), linetype=1, color="white")
psa <- psa + geom_text(aes(x=mstim, y=550, label=paste("mean=", formatC(mstim, 3), sep="")), color="white")
psa <- psa + geom_text(aes(x=mstim, y=150, label="95% HDI"), color="white") + labs(title="PSA enhancer, DHT effect", x="log ff/rn", y="probability mass")
pdf("psaPlusMinusDHT.pdf")
psa
dev.off()


#################################
# SECTION: ANALYSIS of ENHANCER #
#################################

enhsam <- foo$enhEff
n <- NROW(enhsam)
dim(enhsam) <- c(NROW(enhsam)*NCOL(enhsam),1);
enhsam <- data.frame(enheff = enhsam , enh = as.vector(sapply( levels( a$enhancer ) , function(x) rep(x , n))) )
pdf("enhancer_effects.pdf", w=12, h=6)
qplot(enh , enheff , data = enhsam , geom="boxplot") + theme( axis.text.x = element_text( angle = 90 , hjust = 1)) +  labs( title="Basal Enhancer Activity (unstimulated)" , x = "" , y = "log firefly/renilla" )
dev.off()


###################################
# SECTION: ANALYSIS of DHT EFFECT #
###################################

dhtsam <- foo$dhtEff
n <- NROW(dhtsam)
dim(dhtsam) <- c(NROW(dhtsam)*NCOL(dhtsam),1);
dhtsam <- data.frame(dhteff = dhtsam , dht = as.vector(sapply( levels( a$enhancer ) , function(x) rep(x , n))) )

pdf("dht_effects.pdf", w=12, h=6)

ggplot(dhtsam) + geom_boxplot(aes(dht, dhteff)) + theme( axis.text.x = element_text( angle = 90 , hjust = 1)) + ylim( -1 , 4) + labs( title="DHT Effect By Enhancer" , x = "" , y = "log fold change" )

dev.off()


###############################
# SECTION: ANALYSIS of BATCH  #
###############################


batchsam <- subset(foo$batchEff, select=1:6)
n <- NROW( batchsam )
dim( batchsam ) <- c( NROW( batchsam ) * NCOL( batchsam ) , 1 ) ;
batchsam <- data.frame( batcheff=batchsam , batch=c( sapply( levels( as.factor( fflist$batch ) ) , function(x) rep( x , n )) , n )[1:60000])
pdf( "batch_effects.pdf" )
qplot( batch , batcheff , data=batchsam , geom="boxplot" ) + theme( axis.text.x = element_text( angle=90 , hjust = 1 )) + labs( title = "Batch Effect" , x = "" , y = "log(firefly/renilla)")
dev.off()

###############################
PLASMID PREPARATION AND TRANSFECTION EFFECTS

pdf("plasmideff.pdf")
hdi <- HDIofMCMC(foo$sdplas)
ggplot(data.frame(sdsams=foo$sdplas)) + geom_histogram(aes(sdsams)) + labs(title="plasmid effect st dev", x="log fold change", y="probability mass") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
dev.off()

pdf("transfectioneff.pdf")
hdi <- HDIofMCMC(foo$trsd)
ggplot(data.frame(sdsams=foo$trsd)) + geom_histogram(aes(sdsams)) + labs(title="transfection effect st dev", x="log fold change", y="probability mass") + geom_segment(aes(x=hdi[1], xend=hdi[1], y=0, yend=500), linetype=3 ) + geom_segment(aes(x=hdi[2], xend=hdi[2], y=0, yend=500), linetype=3) + geom_text(aes(x=c(hdi), y=rep(550,2), label=formatC(hdi, digits=3)))
dev.off()

