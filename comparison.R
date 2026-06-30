library(hyperinf)
library(hyperdags)
library(ggpubr)
library(ape)

sf = 2

###### illustrative Mk/HMM plot, linear pathway

set.seed(12)
L = 4
n = 32
#n = 16
sim.tree = hyperdags::simulate_accumulation(n, L, dynamics="linear")
data.0p = do.call(rbind, sim.tree$x)[1:n,]
tree.0p = sim.tree$my.tree
rownames(data.0p) = tree.0p$tip.label

plot_hyperinf_data(data.0p, tree.0p)
fit.0pa = hyperinf(data.0p, tree.0p, boot.parallel = 10)
fit.0pb.raw = hyperinf(data.0p, tree.0p, reversible=TRUE)
plot_hyperinf(fit.0pb.raw)
plot_hyperinf_ordering_matrices(list(fit.0pa, fit.0pb.raw))
fit.0pb = bootstrap_mk_fit(fit.0pb.raw)
plot_hyperinf_bubbles(fit.0pb)
#fit.0pc = hyperinf(data.0p, tree.0p, method="hypertraps")

co.plot.2 = plot_hyperinf_compare_orderings(fit.0pa, fit.0pb, p.scale=0.5, thetastep=3, expt.names=c("HMM", "Mk"))

#plot_hyperinf_bubbles(list(fit.0pa, fit.0pb, fit.0pc), p.scale=0.5)


plot.om.abs = plot_hyperinf_ordering_matrices(c(fit.0pa$boots, fit.0pb$boots),
                                              expt.names=c(rep("HMM",11), rep("Mk", 10)), type="absolute")
plot.om.rel = plot_hyperinf_ordering_matrices(c(fit.0pa$boots, fit.0pb$boots),
                                              expt.names=c(rep("HMM",11), rep("Mk", 10)), type="transitions")

co.net.plot = plot_hyperinf_comparative(c(fit.0pa$boots[1:5], fit.0pb$boots[1:5]), style="full",
                                        expt.names=c(rep("HMM", length(fit.0pa$boots[1:5])),
                                                     rep("Mk", length(fit.0pb$boots[1:5]))), threshold = 0.02, bend=2) 

part.1.plot = ggarrange(
  plot_hyperinf_data(data.0p, tree.0p), co.net.plot,
  plot.om.abs, plot.om.rel, labels=c("A", "B", "C", "D"))

png("part-1-plot.png", width=600*sf, height=480*sf, res=72*sf)
print(part.1.plot)
dev.off()

######

set.seed(1)
sim.tree.2 = hyperdags::simulate_accumulation(n, L, dynamics="poisson")
data.0 = do.call(rbind, sim.tree.2$x)[1:n,]
tree.0 = sim.tree.2$my.tree
rownames(data.0) = tree.0$tip.label

plot_hyperinf_data(data.0, tree.0)
fit.0a = hyperinf(data.0, tree.0, boot.parallel = 10)
fit.0b.raw = hyperinf(data.0, tree.0, reversible=TRUE)
fit.0b = bootstrap_mk_fit(fit.0b.raw)
#fit.0c = hyperinf(data.0, tree.0, method="hypertraps")

co.plot.1 = plot_hyperinf_compare_orderings(fit.0a, fit.0b, p.scale=0.5, thetastep=3, expt.names=c("HMM", "Mk"))

plot.om.abs.1 = plot_hyperinf_ordering_matrices(c(fit.0a$boots, fit.0b$boots),
                                              expt.names=c(rep("HMM",11), rep("Mk", 10)), type="absolute")
plot.om.rel.1 = plot_hyperinf_ordering_matrices(c(fit.0a$boots, fit.0b$boots),
                                              expt.names=c(rep("HMM",11), rep("Mk", 10)), type="transitions")

#plot_hyperinf_bubbles(list(fit.0a, fit.0b, fit.0c), p.scale=0.5)

part.2.plot = ggarrange(plot_hyperinf_data(data.0, tree.0), 
          plot_hyperinf_data(data.0p, tree.0p),
          plot.om.abs.1, plot.om.abs,
          plot.om.rel.1, plot.om.rel,
          co.plot.1, co.plot.2, labels=c("Ai","Bi", 
                                         "ii", "ii",
                                         "iii", "iii",
                                         "iv", "iv"),
          ncol = 2, nrow=4)

png("part-2-plot.png", width=600*sf, height=800*sf, res=72*sf)
print(part.2.plot)
dev.off()

######## cross-sectional example

set.seed(1)
n = 4
data.a = matrix(rep(c(0,0,0,0,1, 
                      0,0,0,1,1, 
                      0,0,1,1,1,
                      0,1,1,1,1), n/4), byrow = TRUE, ncol=5, nrow=n)
fit.1a = hyperinf(data.a, boot.parallel = 30)
fit.2a = hyperinf(1-data.a, boot.parallel = 30)
plot.co.a = plot_hyperinf_compare_orderings(fit.1a, fit.2a, thetastep=3, p.scale=0.3)

fig.a = ggarrange(ggarrange(plot_hyperinf_data(data.a),
                            plot_hyperinf_data(1-data.a), nrow=1, labels=c("", "ii")), 
                  plot.co.a, nrow=2, labels=c("", "iii"))

n = 12
data.b = matrix(rep(c(0,0,0,0,1, 
                      0,0,0,1,1, 
                      0,0,1,1,1,
                      0,1,1,1,1), n/4), byrow = TRUE, ncol=5, nrow=n)
fit.1b = hyperinf(data.b, boot.parallel = 30)
fit.2b = hyperinf(1-data.b, boot.parallel = 30)
plot.co.b = plot_hyperinf_compare_orderings(fit.1b, fit.2b, thetastep=3, p.scale=0.3)

fig.b = ggarrange(ggarrange(plot_hyperinf_data(data.b),
                            plot_hyperinf_data(1-data.b), nrow=1, labels=c("", "ii")), 
                  plot.co.b, nrow=2, labels=c("", "iii"))


part.3.plot = ggarrange(fig.a, fig.b, labels=c("Ai", "Bi"))

png("part-3-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(part.3.plot)
dev.off()

ggarrange(plot_hyperinf(fit.1b),
          plot_hyperinf(fit.2b))

####### KpAMR example

# pull from Kp paper
name <- load("~/Dropbox/klebevo/kp-evolution-inference-curate-wip/kleborate-analysis/all_models.Rdata")
country.list <- get(name)

# here we see differences
fit.x = multiple_fits_to_booted_fit(country.list$Gambia)
fit.y = multiple_fits_to_booted_fit(country.list$South_Korea)

kp.amr.plot = plot_hyperinf_compare_orderings(fit.x, fit.y, 
                                         sqrt.trans=TRUE, thetastep=3,
                                         expt.names=c("Gambia", "South Korea"), 
                                         feature.names = fit.x$boots[[1]]$featurenames)

png("part-kp-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(kp.amr.plot)
dev.off()

####### MRO example (not used in current draft)

api.tree = read.tree("mro-ncbi-tree-2025-apicomplexans.nwk")
cil.tree = read.tree("mro-ncbi-tree-2025-ciliophora.nwk")
api.df = read.csv("mro-barcodes-2025-apicomplexans.csv")
cil.df = read.csv("mro-barcodes-2025-ciliophora-1.csv")
api.tree$tip.label = gsub("_", " ", api.tree$tip.label)
cil.tree$tip.label = gsub("_", " ", cil.tree$tip.label)

plot_hyperinf_data(api.df, api.tree)
plot_hyperinf_data(cil.df, cil.tree)

api.fit = hyperinf(api.df, api.tree, losses = TRUE, boot.parallel = 10)
cil.fit = hyperinf(cil.df, cil.tree, losses = TRUE, boot.parallel = 10)

mro.plot = plot_hyperinf_compare_orderings(api.fit, cil.fit,
                                           thetastep=3, p.scale=0.5,
                                           expt.names = c("Apicomplexans", "Ciliates"))

ggarrange(kp.amr.plot, mro.plot, labels=c("A", "B"), widths=c(1.6,1))

