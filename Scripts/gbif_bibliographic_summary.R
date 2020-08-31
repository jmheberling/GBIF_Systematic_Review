##GBIF systematic review of literature 2020
#Bibliographic summaries
#This script contains annotated code to recreate results behind Figures 1,2,4 and supplemental figures 
#see separate script for topic model results (Figure 3 and related results) - "scripts/gbif_literature_topicModel.R"
#last edit 28 Aug 2020

#read in packages--------------------------------------------
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(RColorBrewer) 

#read in bibliographic dataset (GBIF-mediated studies compiled through GBIF literature tracking program)---------------------------
#see Methods and here for description of literature tracking: https://www.gbif.org/literature-tracking
#literature dataset is continuously updated, available here: https://www.gbif.org/resource/search?contentType=literature
#see "Data/GBIF_bibliographic_dataset_metadata.csv" for descriptions of columns
bib<-read_csv("Data/GBIF_bibliographic_dataset.csv") #latest version 16 June 2020 (all years) with data DOI and taxonomic fields
dim(bib)[1] #4153 papers (all years) 
length(bib$year[bib$year>2015]) #2496 papers (recent 2016-19 years only)

#some housekeeping before summarizing data
bib$source<-as.factor(bib$source) #coerce to make journal factor data type
bib$year<-as.integer(bib$year) #make year integer data type
bib$taxonomic.group<-as.factor(bib$taxonomic.group) #make taxonomic groupings factor data type
bib$`tracking source`<- as.factor(bib$`tracking source`) #make tracking source factor data type

#add binary column to denote which papers cite the doi
bib$cite.doi<-rep(0, dim(bib)[1])
bib$cite.doi[bib$`Number of downloads`>0]<-1
sum(bib$cite.doi)#520 papers with at least one data doi cited

#create new dataframe, subsetting to recent studies only (2016-2019 inclusive)
bib.recent<-bib[bib$year>2015,] 

#add binary column for whether studies are "global" or "regional" (i.e., study about >1 GBIF region)
bib.recent$multipleRegions<-rep(0, dim(bib.recent)[1])
bib.recent$multipleRegions[bib.recent$countriesOfCoverage=="None (global or NA)"] <-1
sum(bib.recent$multipleRegions) #1723 studies with no country listed
length(bib.recent$multipleRegions[bib.recent$multipleRegions==0]) #773 single region studies

#summary of data tracking source (for PRISMA diagram, supplemental figure)------
sort(summary(bib$`tracking source`[bib$year>2015])) #only tracked in 2016 onward
dim(bib)[1] #4153 papers (all years); 4035 papers with abstracts and included in topic model
length(bib$year[bib$year>2015]) #2496 papers (recent 2016-19 years only)


#summary of journals published (supplemental figure)----------------------

#Journal frequencies (all years)
JournalFreq<-ddply(bib,c("source"),summarise, no=length(source),no.open=length(source[openAccess==T]))
sum(JournalFreq$no.open)/sum(JournalFreq$no) #30% articles open access (at time of publication)
JournalFreq<-data.frame(JournalFreq) #make dataframe
JournalFreq<-JournalFreq[order(-JournalFreq$no),] #descending order by frequency for barplot aesthetics

#Journal frequencies (recent years only, 2016-19)
JournalFreq.recent<-ddply(bib.recent,c("source"),summarise, no=length(source),no.open=length(source[openAccess==T]))
JournalFreq.recent<-data.frame(JournalFreq.recent) #make dataframe
JournalFreq.recent<-JournalFreq.recent[order(-JournalFreq.recent$no),] #descending order by #papers in journal
sum(JournalFreq.recent$no.open)/sum(JournalFreq.recent$no) #38% open access


#summary of who and where (geography) (Fig. 2)-------------------------------

#Countries of coverage in studies**
#separate country codes into separate rows for country counts (separated by commas in "countriesOfCoverage" for each paper)
countriesCovered<-separate_rows(bib.recent,countriesOfCoverage,sep=",") #recent years only (2016-2019)
#summarize counts by country
countriesCount<-ddply(countriesCovered,.(countriesOfCoverage),summarize,n.coverage=length(countriesOfCoverage))

#GBIF node regions follow Brooks et al. (2016) Scientific Data 3:160007 (see their Fig. 1B "https://www.nature.com/articles/sdata20167/figures/1")
#gbifapi package for extracting GBIF regions from 
#see https://github.com/jhnwllr/gbif_regional_statistics
#install packages built by John Waller
#devtools::install_github("jhnwllr/gbifapi") 
#devtools::install_github("jhnwllr/gbifregionalstats", subdir="gbif_regional_statistics")
library(gbifapi) #load package

#dataframe of GBIF region assignment by countries; see "Data/GBIF_countrycodes.csv" (same as "GBIF_countrycodes" below)
GBIF_countrycodes<-get_gbif_countries() 

#merge list of countries covered (unique country by study) with GBIF region assignment 
GBIF_region_pubs<-merge(countriesCovered[,13],GBIF_countrycodes[,c(2,3)],by.x="countriesOfCoverage",by.y="iso2",all.x=T) 

#summarize unique country coverage by study for each GBIF region
GBIFregion.counts<-ddply(GBIF_region_pubs,.(gbifRegion),summarize,n=length(gbifRegion))
GBIFregion.counts$gbifRegion[8]<-"multiRegion" #these are multiregion studies with no countries listed

#add columns calculating relevant proportions 
GBIFregion.counts$prop.all<-GBIFregion.counts$n/2496*100 #proportion of of recent (2016-19) studies by region (including global studies)
GBIFregion.counts$prop.singleRegion<-GBIFregion.counts$n/(2496-1723)*100 #proportion of studies by region (single region studies only)
GBIFregion.counts$prop.singleRegion[is.na(GBIFregion.counts$gbifRegion)]<-NA #make multi region studies "NA"

GBIFregion.counts #view regional coverage of studies and relevant proportions

#Country-level author affiliations by region**
#similar summary as above, but for country-level affiliation of study authors ("countriesOfResearcher")
#note that below summarizes unique affiations within study, counting number of publications with at least 1 author affiliated with given country (unique author affiliation by study combinations, not distinguishing authors across papers)
length(bib.recent$countriesOfResearcher[!is.na(bib.recent$countriesOfResearcher)]) #all 2496 recent studies are assigned author affiliation(s)

#separate country codes into separate rows for country counts (separated by commas in "countriesOfResearcher" for each paper)
countriesAffil<-separate_rows(bib.recent,countriesOfResearcher,sep=",")
#summarize unique country affiliation counts by country
researcherCount<-ddply(countriesAffil,.(countriesOfResearcher),summarize,n.affiliation=length(countriesOfResearcher)) 

length(countriesAffil$countriesOfResearcher)#4933 unique author affiliation by study combinations (including regional and global studies)
length(countriesAffil$countriesOfResearcher[countriesAffil$multipleRegions==0])#1170 author affiliation by study combinations (single region studies only)

#merge list of country-level authorship (unique country affiliation by study) with GBIF region assignment 
GBIF_region_pubs.affil<-merge(countriesAffil[,12],GBIF_countrycodes[,c(2,3)],by.x="countriesOfResearcher",by.y="iso2",all.x=T)

#summarize unique country-level authorship for each GBIF region
GBIFregion.counts.affil<-ddply(GBIF_region_pubs.affil,.(gbifRegion),summarize,n=length(gbifRegion))

#add row for Antarctica to match study coverage (no authors affiliated with Antarctica)
GBIFregion.counts.affil<-rbind(GBIFregion.counts.affil,c("ANTARCTICA",0))

#make sure region counts are formatted as numbers
GBIFregion.counts.affil$n<-as.numeric(GBIFregion.counts.affil$n) 

#proportions of recent (2016-19) authorship by region (out of all unique author affiliation by study combinations)
GBIFregion.counts.affil$prop.all<-GBIFregion.counts.affil$n/dim(countriesAffil)[1]*100
GBIFregion.counts.affil

#paneled figure of region studied and affiliation by regions (Fig. 2B and Fig. 2C)
par(mfrow=c(2,1),mai=c(.5,0.5,0,0),oma=c(4,4,1,1),xpd=F)
Fig2B<-barplot(sort(GBIFregion.counts$n[-8],decreasing=F),col="white",horiz=F,las=2,ylim=c(0,300),cex.names=1.5)
mtext("# studies using data from:",side = 2, line =4.5,cex=1)
mtext(expression("("*italic("single region studies only")*")"),side = 2, line =3,cex=0.8)
Fig2C<-barplot(GBIFregion.counts.affil$n[c(7,6,1,5,3,2,4)],col="white",horiz=F,las=2,ylim=c(0,2000),cex.names=1.5) #,names.arg=c("Antarctica","Oceania","Africa","N. America","Europe","Asia","Latin America")
mtext("# authors from:",side = 2, line =4.5,cex=1)
mtext(expression("("*italic("unique affiliation by study counts")*")"),side = 2, line =3,cex=0.8)
text(x = Fig2B[,1],y = par("usr")[3] - 100,labels = c("Antarctica","Oceania","Africa","N. America","Europe","Asia","Latin America"),xpd = NA,srt = 45,adj = 0.965,cex = 1)

#To recreate Fig. 2A (world map superimposing country-level authorship vs. coverage) see "Output/Fig2A_GBIF Authorship_by Coverage_Map_Kepler.json"
#Map created using Kepler.gl (open source); go to https://kepler.gl/demo and upload "Fig2A_GBIF Authorship_by Coverage_Map_Kepler.json"

#author affiliation alternative versions (supplemental figures for comparison to Fig. 2) -----
#similar to authorship summary above, but now separating global studies (ie. >1 region covered) and single region studies

#single region studies version**
#separate country codes into separate rows for country counts (separated by commas in "countriesOfResearcher" for each paper)
countriesAffilsingle<-separate_rows(bib.recent[bib.recent$multipleRegions==0,],countriesOfResearcher,sep=",")
ddply(countriesAffilsingle,.(countriesOfResearcher),summarize,n=length(countriesOfResearcher)) 
length(countriesAffilsingle$countriesOfResearcher)#1170 author affiliation by study (single region studies only)
GBIF_region_pubs.affilsingle<-merge(countriesAffilsingle[,12],GBIF_countrycodes[,c(2,3)],by.x="countriesOfResearcher",by.y="iso2",all.x=T)
GBIFregion.counts.affilsingle<-ddply(GBIF_region_pubs.affilsingle,.(gbifRegion),summarize,n=length(gbifRegion))
GBIFregion.counts.affilsingle<-rbind(GBIFregion.counts.affilsingle,c("ANTARCTICA",0))
GBIFregion.counts.affilsingle$n<-as.numeric(GBIFregion.counts.affilsingle$n) #make sure region counts are formatted as numbers

#proportions of recent (2016-19) authorship by region (out of all unique author-publication combinations)
GBIFregion.counts.affilsingle$prop.all<-GBIFregion.counts.affilsingle$n/dim(countriesAffilsingle)[1]*100
GBIFregion.counts.affilsingle

#barplot of authorship counts by region (*single region studies only)
par(mfrow=c(1,1),mai=c(.5,0.5,0,0),oma=c(7,4,1,1),xpd=F)
barplot(GBIFregion.counts.affilsingle$n[c(7,6,1,5,3,2,4)],col="white",horiz=F,las=2,ylim=c(0,2000),names.arg=c("Antarctica","Oceania","Africa","N. America","Europe","Asia","Latin America"),cex.names=1.5)
mtext("# authors from each region",side = 2, line =4.5,cex=1.3)
mtext(expression("("*italic("unique author affiliation by study counts")*")"),side = 2, line =3,cex=0.8)

#global studies version**
#separate country codes into separate rows for country counts (separated by commas in "countriesOfResearcher" for each paper)
countriesAffilglobal<-separate_rows(bib.recent[bib.recent$multipleRegions==1,],countriesOfResearcher,sep=",")
ddply(countriesAffilglobal,.(countriesOfResearcher),summarize,n=length(countriesOfResearcher)) 
length(countriesAffilglobal$countriesOfResearcher)#3763 author affiliation by study (global region studies only)
GBIF_region_pubs.affilglobal<-merge(countriesAffilglobal[,12],GBIF_countrycodes[,c(2,3)],by.x="countriesOfResearcher",by.y="iso2",all.x=T)
GBIFregion.counts.affilglobal<-ddply(GBIF_region_pubs.affilglobal,.(gbifRegion),summarize,n=length(gbifRegion))
GBIFregion.counts.affilglobal<-rbind(GBIFregion.counts.affilglobal,c("ANTARCTICA",0))
GBIFregion.counts.affilglobal$n<-as.numeric(GBIFregion.counts.affilglobal$n) #make sure region counts are formatted as numbers

#proportions of recent (2016-19) authorship by region (out of all unique author-publication combinations)
GBIFregion.counts.affilglobal$prop.all<-GBIFregion.counts.affilglobal$n/dim(countriesAffilglobal)[1]*100
GBIFregion.counts.affilglobal

#barplot of authorship counts by region (*global studies only)
par(mfrow=c(1,1),mai=c(.5,0.5,0,0),oma=c(7,4,1,1),xpd=F)
barplot(GBIFregion.counts.affilglobal$n[c(7,6,1,5,3,2,4)],col="white",horiz=F,las=2,ylim=c(0,2000),names.arg=c("Antarctica","Oceania","Africa","N. America","Europe","Asia","Latin America"),cex.names=1.5)
mtext("# authors from each region",side = 2, line =4.5,cex=1.3)
mtext(expression("("*italic("unique author affiliation by study counts")*")"),side = 2, line =3,cex=0.8)

#summary of GBIF data availability over time, datasets, data DOI citations, etc (Fig. 1)------------------

summary(bib$`Number of downloads`) #number of separate downloads per paper that cited data doi(s)
summary(bib$`Total records in downloads`) #number of records per paper that cited data doi(s)

#spreadsheet of citation counts by dataset
datasets.cited<-read_csv("Data/gbif_dataset_citations_21July2020.csv")
length(datasets.cited$datasetkey[datasets.cited$citations>0]) #26046 datasets cited 
length(datasets.cited$datasetkey[datasets.cited$citations>0])/length(datasets.cited$datasetkey) #48% of datasets with at least one citation

# GBIF datasets sizes by date snapshots, including basis of record, citizen science tag
GBIF_datasets<-read_csv("Data/gbif_dataset_basisofrecord_occurrence_counts_22July2020.csv")
GBIF_datasets$snapshot_date<-as.Date(GBIF_datasets$snapshot_date,format="%m/%d/%y")

#group basis of record into observations vs. specimens (sensu Troudet et al. (2018) Systematic Biology 67:1110-1119)
GBIF_datasets$basis2<-rep(NA,length(GBIF_datasets$occurrence_count))
GBIF_datasets$basis2[GBIF_datasets$basis_of_record=="FOSSIL_SPECIMEN"|GBIF_datasets$basis_of_record=="LIVING_SPECIMEN"|GBIF_datasets$basis_of_record=="MATERIAL_SAMPLE"|GBIF_datasets$basis_of_record=="PRESERVED_SPECIMEN"] <-"Specimen-based"
GBIF_datasets$basis2[GBIF_datasets$basis_of_record=="HUMAN_OBSERVATION"|GBIF_datasets$basis_of_record=="LITERATURE"|GBIF_datasets$basis_of_record=="MACHINE_OBSERVATION"|GBIF_datasets$basis_of_record=="OBSERVATION"] <-"Observation-based"
GBIF_datasets$basis2[GBIF_datasets$basis_of_record=="UNKNOWN"]<-"Unknown"

#total number occurrences across all datsets by snapshots
GBIF_datasets.tot<- GBIF_datasets %>%
  group_by(snapshot_date) %>%
  summarise(total = sum(occurrence_count))

#calculate percent increase all GBIF records from 2007 to 2020
((1570642438-125604611)/125604611)*100

GBIF_datasets.tot$basis2<-rep("total",length(GBIF_datasets.tot$snapshot_date))
GBIF_datasets.tot<-GBIF_datasets.tot[,c(1,3,2)]

#dataset sized by snapshots and by basis of record, counts and proportions by snapshot
GBIF_datasets.byBasis<- GBIF_datasets %>%
  group_by(snapshot_date,basis2) %>%
  summarise(total = sum(occurrence_count)) %>%
  mutate(prop=total/sum(total))

GBIF_datasets.byBasis[GBIF_datasets.byBasis$snapshot_date==as.Date("2020-07-01"),] #85% of records are observation based, 14% specimen based (as of July 2020)

#218831824 specimen based records
#Arino (2010) Biodiversity Informatics 7:81-92 estimates ~1-2 Bn specimens worldwide
#July 2020 - ~219 million specimens-based records in GBIF 
218831824/1000000000;218831824/2000000000 #about 10-20% in GBIF in 2020

#dataset sized by snapshots and by citizen science tag
GBIF_datasets.CitizenScience<- GBIF_datasets %>%
  group_by(snapshot_date,tagged_citizen_science) %>%
  summarise(total = sum(occurrence_count)) %>%
  mutate(prop=total/sum(total))

#65% of records are "citizen science" (July 2020) (1024161718/(1024161718+546480720)), but only account for less than 2% (1.5%) of datasets

#Subset to include most recent GBIF snapshot only
GBIF_datasets.recent<-GBIF_datasets[GBIF_datasets$snapshot_date=="2020-07-01",]
length(unique(GBIF_datasets.recent$datasetkey)) #31099 recent datasets
length(unique(GBIF_datasets$datasetkey)) #35363 total datasets
length(unique(datasets.cited$datasetkey[datasets.cited$citations>0])) #53815 datasets; 26046 cited

GBIF_datasets.recent$datasetkey<-as.factor(GBIF_datasets.recent$datasetkey)

#note some datasets have more than one basis of record represented, so this number doesn't match directly below
GBIF_datasets.recent%>%
  group_by(tagged_citizen_science) %>%
  count(tagged_citizen_science)

length(unique(GBIF_datasets.recent$datasetkey)) #number of unique datasets
length(unique(GBIF_datasets.recent$datasetkey[GBIF_datasets.recent$tagged_citizen_science==TRUE])) #number of datasets tagged as citizen science
length(unique(GBIF_datasets.recent$datasetkey[GBIF_datasets.recent$tagged_citizen_science==FALSE])) #number of datasets tagged as citizen science

#474/31099 datasets are "citizen science"; <2% of datasets
length(unique(GBIF_datasets.recent$datasetkey[GBIF_datasets.recent$tagged_citizen_science==TRUE]))/length(unique(GBIF_datasets.recent$datasetkey))

#combine recent dataset metadata with citation counts
GBIF_datasets.combined<-merge(GBIF_datasets.recent,datasets.cited,by="datasetkey",all.x=T) #note that 1538 datasets have records with more than one basis of record

#totals by dataset, with citations and citizen science tag; calculate citations per dataset size metric
GBIF_datasets.sum<- GBIF_datasets.combined %>%
  group_by(datasetkey,snapshot_date,tagged_citizen_science,citations) %>%
  summarise(total.occurrences = sum(occurrence_count),citationsperrecord=citations/sum(occurrence_count))

GBIF_datasets.sum$datasetkey<-as.factor(GBIF_datasets.sum$datasetkey)

boxplot(GBIF_datasets.sum$total.occurrences) #note ebird is a huge outlier, and long tail 

boxplot(GBIF_datasets.sum$citations) #long tail

#nonparametric tests to avoid issues with nonnormal data (skewed data with key outlier datasets)
#testing differences in number citations for citizen science datasets vs. not citizen science datasets
wilcox.test(citations~tagged_citizen_science,data=GBIF_datasets.sum) #W=6211175, p<0.001
#testing for number citations accounting for dataset size
wilcox.test(citationsperrecord~tagged_citizen_science,data=GBIF_datasets.sum) #W=11514982, p<0.001

#correlation between dataset size and number of ciations
cor.test(GBIF_datasets.sum$total.occurrences,GBIF_datasets.sum$citations,method="spearman") #rho=0.39; p<0.001, df=(32637-2)


GBIF_datasets.sum %>%
  group_by(tagged_citizen_science) %>%
  summarise(
    count = n(),
    median = median(citations, na.rm = TRUE),
    IQR = IQR(citations, na.rm = TRUE)
  )

GBIF_datasets.sum %>%
  group_by(tagged_citizen_science) %>%
  summarise(
    count = n(),
    median = median(citationsperrecord, na.rm = TRUE),
    IQR = IQR(citations, na.rm = TRUE)
  )

#totals by dataset, with citations and citizen science tag, similar to above but add basis of record
GBIF_datasets.sumBasis<- GBIF_datasets.combined %>%
  group_by(datasetkey,snapshot_date,basis2,tagged_citizen_science,citations) %>%
  summarise(total.occurrences = sum(occurrence_count),citationsperrecord=citations/sum(occurrence_count))

GBIF_datasets.sumBasis$datasetkey<-as.factor(GBIF_datasets.sumBasis$datasetkey)
GBIF_datasets.sumBasis$basis2<-as.factor(GBIF_datasets.sumBasis$basis2)

#test between number of citations for specimen- versus observation-based datasets
wilcox.test(citations~basis2,data=GBIF_datasets.sumBasis[GBIF_datasets.sumBasis$basis2!="Unknown",])
wilcox.test(citationsperrecord~basis2,data=GBIF_datasets.sumBasis[GBIF_datasets.sumBasis$basis2!="Unknown",]) #seems when accounting for dataset size, citizen science datasets cited less

GBIF_datasets.sumBasis %>%
  group_by(basis2) %>%
  summarise(
    count = n(),
    median.perrecord = median(citationsperrecord, na.rm = TRUE),
    IQR.perrecord = IQR(citationsperrecord, na.rm = TRUE),
    median.citation = median(citations, na.rm = TRUE),
    IQR.citation = IQR(citations, na.rm = TRUE)
  )

#proportion of recent studies by taxonomic group studied (Fig. 1 inset)
#note this field is only complete for recent studies
tax.proportions<-group_by(bib[bib$year>2015,],taxonomic.group) %>% 
  summarize(n = length(taxonomic.group)) %>% # count per taxonomic group
  mutate(proportion = n / sum(!is.na(bib$taxonomic.group))) 

#Pie charts used in Fig. 1B inset - percentage of studies by organism group studied
colors <- brewer.pal(6,"Set2")
pie(tax.proportions$proportion,labels=c("Fungi (2%)","Invertebrates (17%)","Many (8%)","Other (1%)","Plants (44%)","Vertebrates (27%)"),col=colors,cex=2)#,col=rgb(27,158,119,maxColorValue=255)
symbols(0, 0, circles = 0.4, add=TRUE, bg="white",inches=F) #add white cirlce to make donut chart

#number of papers published per year (and break up those that cite the doi)
no.papers<-ddply(bib,.(year),summarize,papers=length(year),papers.doi=sum(cite.doi))

#plot of papers published over time (GBIF-mediated data use; Fig. 1B) 
par(mfrow=c(1,1),mai=c(.5,0.5,0,0.5),oma=c(3,3,1,4),xpd=F)
plot(no.papers$year,no.papers$papers,type="n",xlab=NA,ylab=NA,xaxt="n",yaxt="n",cex.axis=1.2,xlim=c(2001,2019))
lines(c(2001,2002,no.papers$year),c(0,0,no.papers$papers),type="l",lwd=6,col="darkgrey")
lines(no.papers$year[no.papers$year>2015],no.papers$papers[no.papers$year>2015],type="l",lwd=6,col="black")
lines(c(2016,no.papers$year[no.papers$year>2015]),c(0,no.papers$papers.doi[no.papers$year>2015]),type="l",lwd=3,lty=2,col="black")
axis(side=1,at=c(2001:2019),las=2,cex.axis=1.2) #could tilt year labels with text(,srt=45)
axis(side=2,at=seq(from=0,to=800,by=100),las=2,cex.axis=1.2) #could tilt year labels with text(,srt=45)

mtext("Year published",side = 1, line =4.2,cex=1.5)
mtext("# GBIF-mediated publications",side = 2, line =3.2,cex=1.5)

#number of records available via GBIF over time (GBIF-mediated data growth; Fig. 1A)

#spreadsheet of occurrence records available via GBIF over time (data snapshots)
#counts by basis of record
GBIF_occ<-read_csv("Data/occ_kingdom_basisOfRecord.csv")
GBIF_occ$snapshot<-as.Date(GBIF_occ$snapshot,format="%m/%d/%y") #format date

#group basis of record into observations vs. specimens (sensu Troudet et al. (2018) Systematic Biology 67:1110-1119)
GBIF_occ$basis2<-rep(NA,length(GBIF_occ$occurrenceCount))
GBIF_occ$basis2[GBIF_occ$basisOfRecord=="FOSSIL_SPECIMEN"|GBIF_occ$basisOfRecord=="LIVING_SPECIMEN"|GBIF_occ$basisOfRecord=="MATERIAL_SAMPLE"|GBIF_occ$basisOfRecord=="PRESERVED_SPECIMEN"] <-"Specimen-based"
GBIF_occ$basis2[GBIF_occ$basisOfRecord=="HUMAN_OBSERVATION"|GBIF_occ$basisOfRecord=="LITERATURE"|GBIF_occ$basisOfRecord=="MACHINE_OBSERVATION"|GBIF_occ$basisOfRecord=="OBSERVATION"] <-"Observation-based"
GBIF_occ$basis2[GBIF_occ$basisOfRecord=="UNKNOWN"]<-"Unknown"

#summarize by date (total)
GBIF_occ.tot<- GBIF_occ %>%
  group_by(snapshot) %>%
  summarise(total = sum(occurrenceCount))

GBIF_occ.tot$basis2<-rep("total",length(GBIF_occ.tot$snapshot))
GBIF_occ.tot<-GBIF_occ.tot[,c(1,3,2)] #reorder columns

#summarize by data and basis of record
GBIF_occ.byBasis<- GBIF_occ %>%
  group_by(snapshot,basis2) %>%
  summarise(total = sum(occurrenceCount))

#plot GBIF-mediated data growth over time (specimen-based, observation-based, and total occurrences)
par(mfrow=c(1,1),mai=c(.5,0.5,0,0.5),oma=c(3,4,1,4),xpd=F)
plot(GBIF_occ.tot$snapshot,GBIF_occ.tot$total/1000000,type="l",xlab=NA,ylab=NA,cex.axis=1.2,xaxt="n",las=2,lwd=4)
axis.Date(side=1,at=seq(as.Date("2001",format="%Y"),as.Date("2019",format="%Y"),by="year"),las=2,cex.axis=1.2) 
lines(GBIF_occ.byBasis$snapshot[GBIF_occ.byBasis$basis2=="Observation-based"],GBIF_occ.byBasis$total[GBIF_occ.byBasis$basis2=="Observation-based"] /1000000,lty=5,lwd=4,col="azure4")
lines(GBIF_occ.byBasis$snapshot[GBIF_occ.byBasis$basis2=="Specimen-based"],GBIF_occ.byBasis$total[GBIF_occ.byBasis$basis2=="Specimen-based"] /1000000,lty=3,lwd=4,col="cornsilk4")

mtext("Year",side = 1, line =4.2,cex=1.5)
mtext("# occurrences accessible via GBIF",side = 2, line =4.7,cex=1.5)
mtext(expression("("*italic("in millions")*")"),side = 2, line =3.2,cex=1)

#breakdown of GBIF-mediated data by coarse taxonomic groups (as of 20 July 2020)
#GBIF kingdomKey: 1=Animalia; 2=Archaea; 3=Bacteria; 4=Chromista; 5=Fungi; 6=Plantae; 7=Protozoa, 8=virus, 0 = doubtful
#using gbifapi (already loaded above)
#these numbers change frequently (numbers in comments were used in publication; as of 20 July 2020)
sum(get_occ_count_from_taxonkeys(c(2,3,4,7,8))) #28350481 = "other" -- 
get_occ_count_from_taxonkeys(c(5)) #19585137 = "fungi"
get_occ_count_from_taxonkeys(c(6)) #296101351 = "plant"
get_occ_count_from_taxonkeys(c(44)) #1066704222 = vertebrates (phylum Chordata)
get_occ_count_from_taxonkeys(c(1))-get_occ_count_from_taxonkeys(c(44)) #156478132 = invertebrates  ( NOT phylum Chordata)

total.kingdom = sum(get_occ_count_from_taxonkeys(c(1:8))) #total 1567219352 assigned to a kingdom

gbif.tax<-data.frame("taxon"=c("Fungi","Invertebrates","Other","Plants","Vertebrates"),n = c(19585137,156478132,28350481,296101351,1066704222))
gbif.tax$prop<-gbif.tax$n/1567219352

#Pie charts used in Fig. 1A inset - percentage of studies by organism group studied
par(mfrow=c(1,1),mai=c(.5,0.5,0,0.5),oma=c(3,4,1,4),xpd=F)
colors <- brewer.pal(6,"Set2")
pie(gbif.tax$prop,labels=c("Fungi (1%)","Invertebrates (10%)","Other (2%)","Plants (19%)","Vertebrates (68%)"),col=colors[c(1,2,4,5,6)],cex=2)#,col=rgb(27,158,119,maxColorValue=255)
symbols(0, 0, circles = 0.4, add=TRUE, bg="white",inches=F)

#Map of Science summaries by disciplines (Fig. 4)------------------
#data from Sci2 tool available at Sci2 Team. (2009). Science of Science Tool (Sci2). Indiana University and SciTech Strategies, https://sci2.cns.iu.edu.
#map in manuscript generated using Sci2 tool
#see Borner et al (2012) PLoS ONE 7:e39464 for details, including journal subdiscipline and major discipline assignments in their supplemental information
#see "Data/GBIF_WoS_allYears 4June2020.txt" for exported full Web of Science records (including references) for GBIF-enabled studies indexed in Web of Science
#see "Output/GBIF_map_of_science_journals_mapped.csv" for list of journals mapped, number of studies per journal


