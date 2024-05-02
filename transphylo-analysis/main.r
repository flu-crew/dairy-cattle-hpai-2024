#!/usr/bin/env Rscript

if (!dir.exists("./.RLibs")){
dir.create("./.RLibs", recursive = TRUE)
}

.libPaths("./.RLibs")  # add to the path

packages <- c('ape', 'argparse', 'lubridate', 'TransPhylo', 'coda', 'lattice')
install.packages(setdiff(packages, rownames(installed.packages())), lib="./.RLibs" , repos = "https://cloud.r-project.org")

library(TransPhylo)
library(argparse)
library(lubridate)
library(ape)
library(coda)
library(lattice)
set.seed(22)

parser <- ArgumentParser(description='Run transmission analyses with TransPhylo')
parser$add_argument('treefile', metavar='T', type="character", help='Treefile to analyze')
parser$add_argument('-d', '--date', metavar='date', type="character", help='Date last sampled')
parser$add_argument('--g-mean', metavar='gmean', type="double", help='Gama distribution mean of generation time (in years)')
parser$add_argument('--s-mean', metavar='smean', type="double", help='Gama distribution mean of sampling time density (in years)')
parser$add_argument('--g-std', metavar='gstd', type="double", help='Gama distribution standard deviation of generation time (in years)')
parser$add_argument('--s-std', metavar='sstd', type="double", help='Gama distribution standard deviation of sampling time density(in years)')
parser$add_argument('-i', '--iter', metavar='iter', type="integer", help='Number of MCMC iterations')
args <- parser$parse_args()

# Gamma distribution shape and scale
g_gamma_shape<-args$g_mean/(args$g_std^2)
g_gamma_scale<-args$g_mean/g_gamma_shape

# Get date of most recently sampled taxa
date<-decimal_date(ymd(args$date))

# Read nexus file 
phy<-read.nexus(args$treefile)

# Binarize tree
phy <- multi2di(phy, tol=1e-08) # remove multifurcations
phy$edge.length <- pmax(phy$edge.length,1/365)

# Set window height based on num taxa. uses text size 12 
device_height<-((12/72)*Ntip(phy))+4

# plot tree
ptree<-ptreeFromPhylo(phy,dateLastSample=date)
png("phylotree.png", height=device_height*2, width=15, units="in", res=62)
plot(ptree, x.lim=15)
dev.off()

res<-inferTTree(ptree, mcmcIterations=args$iter, w.mean=args$g_mean, w.std=args$g_std, ws.mean=args$s_mean, ws.std=args$s_std, dateT=Inf, thinning=10)

# Print the ESS (effective sample size) of the paramters 
mcmc=convertToCoda(res)
print(effectiveSize(mcmc))

# Print mean and 95% credible interval of parameters
print(res)


# Plot medoid transmission tree. See TransPhylo for more details
med=medTTree(res)
png("Coloured_ttree_detailed.png", height=device_height, width=device_height*5, units="in", res=72)
plot(extractTTree(med), type='detailed', w.shape=g_gamma_shape, w.scale=g_gamma_scale)
dev.off()

png("Coloured_ttree.png", height=device_height, width=38, units="in", res=72)
plot(med)
dev.off()

# Probability of pairwise transmission
mat=computeMatWIW(res)
png("pairwise_transmission_probs.png", height=device_height, width=device_height, units="in", res=72)
levelplot(mat,xlab='',ylab='', default.scales=list(x=list(rot=45)))
dev.off()

# Number of intermediates in the transmission chain
mat=computeMatTDist(res)
png("num_transmission_intermediates.png", height=device_height, width=device_height, units="in", res=72)
levelplot(mat,xlab='',ylab='', default.scales=list(x=list(rot=45)))
dev.off()

# Distribution of detected cases, generation times, and sampling time
png("Incident_cases.png", height=10, width=10, units="in", res=72)
a<-getIncidentCases(res, numBins=24, dateT=2024.2, show.plot=T)
dev.off()

png("generation_time_distribution.png", height=10, width=10, units="in", res=72)
getGenerationTimeDist(res,show.plot = T)
dev.off()

png("sampling time distribution.png", height=10, width=10, units="in", res=72)
getSamplingTimeDist(res,show.plot = T)
dev.off()

# Plot MCMC Traces
png("mcmc_traces.png")
plot(res)
dev.off()

