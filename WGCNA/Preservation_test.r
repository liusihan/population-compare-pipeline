#######Calculation of module preservation
load("../ChineseHumanBrain/Chinese-Signed-networkConstruction-auto.RData")
load("../BrainGVEx/Caucasian-Unsigned-networkConstruction-auto.RData")
load("../overlapgenes.rdata")
load("AA-Unsigned-networkConstruction-auto.RData")
color <- moduleColors
#colorGVX <- GVXinformatio$module
setLabels = c("AA","CNS");
multiExpr = list(AA = list(data = t(AA_datExpr[genes,])),CNS = list(data = t(CNS_datExpr[genes,])));
multiColor = list(AA = color)

#ref=AA
system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
nPermutations = 200,
randomSeed = 1,
quickCor = 0,
verbose = 3)
} );
# Save the results
save(mp, file = "AA-UnSigned-ref-CNS-test-modulePreservationin.RData");

load("CAUC-Unsigned-ref-AA-test-modulePreservationin.RData")
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )


# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
pdf(file="CAUCOnly-Unsigned-AA-test-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10) min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.4,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();

# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
pdf(file="CAUCOnly-Unsigned-AA-test-modulePreservation.pdf", wi=10, h=5)
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
min = min(statsZ[plotMods, s], na.rm = TRUE);
max = max(statsZ[plotMods, s], na.rm = TRUE);
if (min > -max/5) min = -max/5
plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
main = colnames(statsZ)[s],
cex = 1.7,
ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
xlim = c(20, 1000))
labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()
