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

co.plot.2 = plot_hyperinf_compare_orderings(fit.0pa, fit.0pb, p.scale=0.5, thetastep=3, expt.names=c("HyperHMM", "HyperMk"))

#plot_hyperinf_bubbles(list(fit.0pa, fit.0pb, fit.0pc), p.scale=0.5)


plot.om.abs = plot_hyperinf_ordering_matrices(c(fit.0pa$boots, fit.0pb$boots),
                                              expt.names=c(rep("HyperHMM",11), rep("HyperMk", 10)), type="absolute")
plot.om.rel = plot_hyperinf_ordering_matrices(c(fit.0pa$boots, fit.0pb$boots),
                                              expt.names=c(rep("HyperHMM",11), rep("HyperMk", 10)), type="relative")

co.net.plot = plot_hyperinf_comparative(c(fit.0pa$boots[1:5], fit.0pb$boots[1:5]), style="full",
                                        expt.names=c(rep("HMM", length(fit.0pa$boots[1:5])),
                                                     rep("Mk", length(fit.0pb$boots[1:5]))), threshold = 0.02, bend=2) 

net.plot.1 = plot_hyperinf_comparative(fit.0pa$boots[1:10], style="limited",
                                       expt.names=rep("HMM", length(fit.0pa$boots[1:10])), 
                                       threshold = 0.02, bend=2, label_size = 2.5) 
net.plot.2 = plot_hyperinf_comparative(fit.0pb$boots[1:10], style="limited",
                                       expt.names=rep("Mk", length(fit.0pa$boots[1:10])), 
                                       threshold = 0.02, bend=2, label_size = 2.5) 

co.net.plot.alt = ggarrange(
net.plot.1 + 
  scale_edge_alpha_continuous(range=c(0.5,1)) + 
  scale_edge_width_continuous(limits = c(0.9,1), range=c(1,2)) +
  guides(edge_colour = "none") +
  scale_edge_color_manual(values = rep("steelblue", 10)) +
  labs(edge_alpha = "Resamples", edge_width = "Flux"),
net.plot.2 + 
  scale_edge_alpha_continuous(range=c(0.1,1)) + 
  scale_edge_width_continuous(limits = c(0,3), range=c(1,2))+
  guides(edge_colour = "none") +
  scale_edge_color_manual(values = rep("tomato", 10)) +
labs(edge_alpha = "Resamples", edge_width = "Flux"),
labels = c("B i", "ii"),
widths=c(1,2)
)

part.1.plot = ggarrange(
  plot_hyperinf_data(data.0p, tree.0p), co.net.plot,
  plot.om.abs, plot.om.rel, labels=c("A", "B", "C", "D"))

recolor = scale_fill_manual(values = c("steelblue", "tomato"))

part.1.plot.alt = ggarrange(
  ggarrange(plot_hyperinf_data(data.0p, tree.0p), co.net.plot.alt, widths=c(0.75,2), labels=c("A", "")),
  ggarrange(plot.om.abs + recolor, plot.om.rel + recolor, labels=c("C", "D")),
  nrow=2)

png("part-1-plot-alt.png", width=700*sf, height=550*sf, res=72*sf)
print(part.1.plot.alt)
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

co.plot.1 = plot_hyperinf_compare_orderings(fit.0a, fit.0b, p.scale=0.5, thetastep=3, expt.names=c("HyperHMM", "HyperMk"))

plot.om.abs.1 = plot_hyperinf_ordering_matrices(c(fit.0a$boots, fit.0b$boots),
                                              expt.names=c(rep("HyperHMM",11), rep("HyperMk", 10)), type="absolute")
plot.om.rel.1 = plot_hyperinf_ordering_matrices(c(fit.0a$boots, fit.0b$boots),
                                              expt.names=c(rep("HyperHMM",11), rep("HyperMk", 10)), type="relative")

#plot_hyperinf_bubbles(list(fit.0a, fit.0b, fit.0c), p.scale=0.5)

part.2.plot = ggarrange(plot_hyperinf_data(data.0, tree.0), 
          plot_hyperinf_data(data.0p, tree.0p),
          plot.om.abs.1 + recolor, plot.om.abs + recolor,
          plot.om.rel.1 + recolor, plot.om.rel + recolor,
          co.plot.1 + recolor, co.plot.2 + recolor, labels=c("Ai","Bi", 
                                         "ii", "ii",
                                         "iii", "iii",
                                         "iv", "iv"),
          heights = c(1,1,1,1.1),
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

recolor2 = scale_fill_manual(values=c("steelblue", "tomato", "white"))
nominorx = scale_x_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)
nominory = scale_y_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)

plot.co.a = plot_hyperinf_compare_orderings(fit.1a, fit.2a, 
                                            expt.names = c("i", "ii"),
                                            thetastep=3, p.scale=0.3)

fig.a = ggarrange(ggarrange(plot_hyperinf_data(data.a),
                            plot_hyperinf_data(1-data.a), nrow=1, labels=c("", "ii")), 
                  plot.co.a + nominorx + nominory +
                    recolor2 + labs(fill = "Dataset"), 
                  nrow=2, labels=c("", "iii"))

n = 12
data.b = matrix(rep(c(0,0,0,0,1, 
                      0,0,0,1,1, 
                      0,0,1,1,1,
                      0,1,1,1,1), n/4), byrow = TRUE, ncol=5, nrow=n)
fit.1b = hyperinf(data.b, boot.parallel = 30)
fit.2b = hyperinf(1-data.b, boot.parallel = 30)
plot.co.b = plot_hyperinf_compare_orderings(fit.1b, fit.2b, 
                                            expt.names = c("i", "ii"),
                                            thetastep=3, p.scale=0.3)

fig.b = ggarrange(ggarrange(plot_hyperinf_data(data.b),
                            plot_hyperinf_data(1-data.b), nrow=1, labels=c("", "ii")), 
                  plot.co.b + nominorx + nominory +
                    recolor2 + labs(fill = "Dataset"), 
                  nrow=2, labels=c("", "iii"))


part.3.plot = ggarrange(fig.a, fig.b, labels=c("Ai", "Bi"))

png("part-3-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(part.3.plot)
dev.off()

ggarrange(plot_hyperinf(fit.1b),
          plot_hyperinf(fit.2b))


