#CLEAN VERSION
# scripted by Brian Leung, Nov 25, 2020 (brian.leung2@mcgill.ca)
#script to process the LPI data, to be consistent with the published analysis
#also, see ids_used_in_analysis.csv for the populations used in the analysis
# please cite:  Leung, B., Hargreaves AL, Greenberg, DA, McGill, B, Dornelas, M, Freeman, R. 2020. Clustered versus catastrophic global vertebrate declines. Nature: https://doi.org/10.1038/s41586-020-2920-6

cyr=3 #the year column
cN=4 #the pop size column
#THESE FILES ARE SAMPLES: THE DATA NEEDS TO BE DOWNLOADED FROM LPI INTO THIS FORM
tot=readRDS("LPI_dat_sample.rds") 
info<-readRDS("info1.rds")


tot=tot[tot[,cyr]>1969 & tot[,cyr]<2015,] 

#there are some with a single value or where all values are zero - we can't use those ones

pos=by(tot,tot$id,function(x){
	if(nrow(x)<2 || all(x[,cN]==0)) #return the indices 
		return(x$id[1])
	return(NA)
})
pos=pos[!is.na(pos)]
tot=tot[!(tot$id %in% pos),]

info=info[,c(2,4,7,8,18,20)]

info$Class=as.character(info$Class) 
info$System=as.character(info$System) 
info$Realm=as.character(info$Realm) 
#combine groups into the 57 combos used in LPI (2018)
#indo-pacific includes Oceania, Australasia, and Indo-Malayan
pos=which(info$Realm %in% c("Oceania","Australasia","Indo-Malayan"))
info$Realm[pos]="indo-Pacific"
info=info[info$Realm != "Antarctic",]
#exclude Plantae
info=info[info$Class !="Plantae",]
info=info[!is.na(info$Class),]
#combine reptiles and amphibians
info$Class[info$Class %in% c("Amphibia","Reptilia")]="Herps"
#then everything else unknown is a fish
info$Class[!info$Class %in% c("Aves","Mammalia","Herps")]="Fish"

info$combo=paste(info$Class,info$System,info$Realm)

#NOTE - SEARCH FOR AND REMOVE DUPLICATES - WHERE BOTH COUNTRY LEVEL AND LOCAL ESTIMATES ARE PROVIDED. WE REMOVED THE HIGHER LEVEL AND KEPT THE SMALLER ONES
#included is a current file with ids for replicates
repl<-read.csv("remv_aggr_pops.csv")
pos=which(info$id %in% repl[,2])
info=info[-pos,]

#force both datafiles to match
info=info[which(info$id %in% unique(tot$id)),]
tot=tot[which(tot$id %in% unique(info$id)),]
growth_rates=by(tot,tot$id,function(x) #this provides growth rates
{
	x=x[order(x$Year),]
	if(length(which(x[,cN]==0))>0) #need to take care of zeros. Change entire time series, adding 1% of mean
	{
		x[,cN]=x[,cN]+0.01*mean(x[,cN])
	}
	dur=x$Year[-1]-x$Year[-nrow(x)] 
	gr<-log(x[-1,cN]/x[-nrow(x),cN])/dur #for non-consecutive years, take mean growth
	return(gr)})

saveRDS(growth_rates,"growth_rates.rds")
saveRDS(info,"info_processed.rds")
