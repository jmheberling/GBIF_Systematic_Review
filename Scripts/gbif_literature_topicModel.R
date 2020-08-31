#GBIF systematic review of literature 2020
#script created 14 Apr 2020; last edit 6 Aug 2020
#structural topic model code using STM package (see Roberts et al. [2014] American Journal of Political Science 58:1262-1272 and Roberts, M.E., Stewart, B.M. & D. Tingley [2019] Journal of Statistical Software 91:1-40 doi: 10.18637/jss.v091.i02)
#code below includes model selection, validation, summaries, and visualizations
#see Methods and Appendix S1 for additional details
#**note that topic numbers as numbered in R output are NOT the same as numbered in paper (topics are numbered by descending frequency in paper for easier reader interpretation)


#read in packages--------------------------------------------
library(stm) # structural topic model package (stm v1.3.5)
library(readr);library(tidyr); library(plyr) #for data handling

#read in bibliographic dataset (GBIF-mediated studies compiled through GBIF literature tracking program)---------------------------
#see Methods and here for description of literature tracking: https://www.gbif.org/literature-tracking
#literature dataset is continuously updated, available here: https://www.gbif.org/resource/search?contentType=literature
#see "Data/GBIF_bibliographic_dataset_metadata.csv" for descriptions of columns
bib<-read_csv("Data/GBIF_bibliographic_dataset.csv") #latest version 16 June 2020 (all years) with data DOI and taxonomic fields
dim(bib)[1] #4153 papers (all years) 
length(bib$year[bib$year>2015]) #2496 papers (recent 2016-19 years only)

#some housekeeping
bib$source<-as.factor(bib$source) #make journal factor
bib$year<-as.integer(bib$year) #make year integer

#create new columns with concatenated fields for analysis
bib$Abstract.Title<-paste(bib$title,bib$abstract,sep=". ") #add title as first sentence of abstract for analysis
bib$Abstract.Title.Keywords<-paste(bib$Abstract.Title,bib$keywords,sep=". ") #also add keywords

#review sample sizes
ddply(bib,.(year),summarize,N=length(year)) #sample size by year; note increase from 2016; few papers in 2003-2007

ddply(bib[!is.na(bib$abstract),],.(year),summarize,N=length(year)) #sample size by year WITH ABSTRACTS only (and therefore included in models)

#Remove records (with reasons)
bib<-bib[bib$id!="0851e8cc-b4d3-3617-961c-36d5f49199e3",] #remove one entry with non-english title AND abstract
bib<-bib[!is.na(bib$abstract),] #Remove records without abstract for topic models (118 papers)

#2016-2019 determined a priori as "recent" based on period when data DOIs were routinely generated with GBIF downloads
dim(bib[bib$year<2016,])[1] #1539 papers in analysis pre-2016
dim(bib[bib$year>2015,])[1] #2496 papers in analysis 2016-2019

dim(bib[bib$year>=2016,])[1]/dim(bib)[1] #2016-19 accounts for 62% of dataset

#Data ready for modeling--------------
#Pipeline for topic structural topic modeling: 
#1.INGEST--> 2.PREPARE --> 3.1EVALUATE & 3.2ESTIMATE --> 4.UNDERSTAND --> 5.VISUALIZE
#following: Roberts, M.E., Stewart, B.M. & D. Tingley [2019] Journal of Statistical Software 91:1-40 doi: 10.18637/jss.v091.i02

#1) INGEST dataset-----------
#stemming, remove white space, etc.: 
#could conisder adding additional stopwords with no useful semantic meaning using customstopwords=c("")
processed <- textProcessor(bib$Abstract.Title.Keywords,metadata=bib) #including abstact, title, and keywords 

#2) PREP dataset--------
#Visualize effect of threshold - that is, number of docs a word must be present in for inclusion (a user-defined parameter!)
#for example, default of lower.thresh=1 means that words which appear in only one document will be dropped
plotRemoved(processed$documents, lower.thresh = seq(1, 500, by = 10)) #visualize numbers of documents, words, and tokens removed at different thresholds

out <- prepDocuments(processed$documents, processed$vocab,processed$meta, lower.thresh = 10) #set lower threshold no. of docs that must contain a word for it to remain in model;upper threshold to remove those words that are too common

#If documents are removed, which ones?
#**Note that if docs are removed, index of docs are different from original dataset (important when interpreting topic meaning and want to read top abstracts associated with each topic)
out$docs.removed #index of removed docs (because they contained no words within threshold defined above)

#3.1) EVALUATE (metrics to assist in model selection)-------
#Compare models with different number of topics, calculate exclusivity and semantic coherence:
#Note this is computationally expensive (takes several hours to run)
#Done in conjunction with ESTIMATE (manually assess model meaning/usefulness with different number of topics modeled)
#Can also compare different model structures (e.g. no covariates) and text inputs (e.g. years, excluding keywords)

#Model structure with linear year covariate (no smoother) (as presented in manuscript)
#here, compare models with 5 to 75 number of topics (K)
Kcompare3<-searchK(out$documents, out$vocab, K = c(5:75), max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))

#model diagnostics across models that vary by K 
plot(Kcompare3)

#plot exclusivity and semantic coherence (see Supplemental Appendix S1)
plot(Kcompare3$results$exclus~Kcompare3$results$semcoh,type="n",ylab="Exclusivity",xlab="Semantic Coherence")
text(Kcompare3$results$exclus~Kcompare3$results$semcoh,label=c(Kcompare3$results$K),col="red")
text(Kcompare3$results$exclus[Kcompare3$results$K==25]~Kcompare3$results$semcoh[Kcompare3$results$K==25],label=c(Kcompare3$results$K[Kcompare3$results$K==25]),col="blue")
#zoom in on topics <40
plot(Kcompare3$results$exclus[Kcompare3$results$K<50]~Kcompare3$results$semcoh[Kcompare3$results$K<50],type="n",ylab="Exclusivity",xlab="Semantic Coherence")
text(Kcompare3$results$exclus[Kcompare3$results$K<50]~Kcompare3$results$semcoh[Kcompare3$results$K<50],label=c(Kcompare3$results$K[Kcompare3$results$K<50]),col="red")

#with linear year covariate (no smoother) - run 2, for comparison with Kcompare3, since not deterministic
Kcompare4<-searchK(out$documents, out$vocab, K = c(5:75), max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
plot(Kcompare4)

#compare two different runs of each model with K ranging from 5-75
plot(Kcompare4$results$exclus~Kcompare4$results$semcoh,type="n",ylab="Exclusivity",xlab="Semantic Coherence")
text(Kcompare4$results$exclus~Kcompare4$results$semcoh,label=c(Kcompare4$results$K),col="red")
text(Kcompare3$results$exclus~Kcompare3$results$semcoh,label=c(Kcompare3$results$K),col="blue") #add prev run

#note that exclusivity and semantic coherence is jointly maximized around 20-30 topics
#Below are subsets of these models for closer inspection and user validation

#3.2) ESTIMATE model (year covariate)--------
#run models with different number of topics; which one is "best" requires user validation, interpretation
#if no covariates included, the model reduces to a (fast) implementation of the correlated topic model (Blei and Lafferty 2007)
#could include year covariate: ", prevalence =~ s(year)," for spline (nonlinear), but in this case (limited years) assume linear
#here (and presented in manuscript), year of publication is included as a covariate
#run many models that vary in number of topics
stm10 <- stm(out$documents, out$vocab, K = 10, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm15 <- stm(out$documents, out$vocab, K = 15,  max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm16 <- stm(out$documents, out$vocab, K = 16,  max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm17 <- stm(out$documents, out$vocab, K = 17,  max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm18 <- stm(out$documents, out$vocab, K = 18,  max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm19 <- stm(out$documents, out$vocab, K = 19,  max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm20 <- stm(out$documents, out$vocab, K = 20, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm21 <- stm(out$documents, out$vocab, K = 21, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm22 <- stm(out$documents, out$vocab, K = 22, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm23 <- stm(out$documents, out$vocab, K = 23, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm24 <- stm(out$documents, out$vocab, K = 24, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm25 <- stm(out$documents, out$vocab, K = 25, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm26 <- stm(out$documents, out$vocab, K = 26, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm27<- stm(out$documents, out$vocab, K = 27, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm28<- stm(out$documents, out$vocab, K = 28, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm29<- stm(out$documents, out$vocab, K = 29, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm30 <- stm(out$documents, out$vocab, K = 30, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm31 <- stm(out$documents, out$vocab, K = 31, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm32 <- stm(out$documents, out$vocab, K = 32, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm33 <- stm(out$documents, out$vocab, K = 33, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm34 <- stm(out$documents, out$vocab, K = 34, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm35 <- stm(out$documents, out$vocab, K = 35, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm40 <- stm(out$documents, out$vocab, K = 40, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm50 <- stm(out$documents, out$vocab, K = 50, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))
stm100 <- stm(out$documents, out$vocab, K = 100, max.em.its = 150, data = out$meta, init.type = "Spectral", prevalence =~ (year))

#4) UNDERSTAND topic meaning---------------------
#below code repeated for each model compared (to assess topic meanings across models) following the same code structure:

#extract top 25 abstracts for each topic; for each model with different number of topics
#add column with top 15 words separated by comma for each topic
#add column with document level theta (proportion of words assigned to topic)
#write to separate files for each model for manual review (interpretation and validation - i.e., do topics make sense and have human interpretive value?)
#carefully read top words and top abstracts for model to validate model (do topics make sense?) and model selection (which model is most appropriate?)

#note that 25 topic model was presented in paper; 20 and 30 topic models were presented in supplemental information for comparison
#**note that topic numbers as numbered in output are NOT the same as they are numbered in paper (topics are numbered by descending frequency in paper for easier reader interpretation)


#10 topics*******
number_of_topics=10
stm10thoughts0<-findThoughts(stm10,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm10thoughts<-findThoughts(stm10,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm10thoughts<-out$meta[data.frame(unlist(stm10thoughts))[,1],]
stm10thoughts$topic <-sort(rep(1:number_of_topics,25))

stm10theta <-data.frame(stm10$theta[data.frame(unlist(stm10thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm10theta)[1]){
  doc.thetas<-stm10theta[i,stm10thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm10thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm10_words<-data.frame((labelTopics(stm10, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm10_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm10thoughts=merge(stm10thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm10thoughts=stm10thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm10thoughts,path="Output/Top_abstracts/stm10_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm10, n = 15)$prob)),path="Output/Top_words/stm10_words.csv")

#15 topics*******
number_of_topics=15
stm15thoughts0<-findThoughts(stm15,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm15thoughts<-findThoughts(stm15,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm15thoughts<-out$meta[data.frame(unlist(stm15thoughts))[,1],]
stm15thoughts$topic <-sort(rep(1:number_of_topics,25))

stm15theta <-data.frame(stm15$theta[data.frame(unlist(stm15thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm15theta)[1]){
  doc.thetas<-stm15theta[i,stm15thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm15thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm15_words<-data.frame((labelTopics(stm15, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm15_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm15thoughts=merge(stm15thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm15thoughts=stm15thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm15thoughts,path="Output/Top_abstracts/stm15_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm15, n = 15)$prob)),path="Output/Top_words/stm15_words.csv")

#16 topics*******
number_of_topics=16
stm16thoughts0<-findThoughts(stm16,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm16thoughts<-findThoughts(stm16,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm16thoughts<-out$meta[data.frame(unlist(stm16thoughts))[,1],]
stm16thoughts$topic <-sort(rep(1:number_of_topics,25))

stm16theta <-data.frame(stm16$theta[data.frame(unlist(stm16thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm16theta)[1]){
  doc.thetas<-stm16theta[i,stm16thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm16thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm16_words<-data.frame((labelTopics(stm16, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm16_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm16thoughts=merge(stm16thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm16thoughts=stm16thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm16thoughts,path="Output/Top_abstracts/stm16_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm16, n = 15)$prob)),path="Output/Top_words/stm16_words.csv")

#17 topics*******
number_of_topics=17
stm17thoughts0<-findThoughts(stm17,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm17thoughts<-findThoughts(stm17,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm17thoughts<-out$meta[data.frame(unlist(stm17thoughts))[,1],]
stm17thoughts$topic <-sort(rep(1:number_of_topics,25))

stm17theta <-data.frame(stm17$theta[data.frame(unlist(stm17thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm17theta)[1]){
  doc.thetas<-stm17theta[i,stm17thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm17thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm17_words<-data.frame((labelTopics(stm17, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm17_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm17thoughts=merge(stm17thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm17thoughts=stm17thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm17thoughts,path="Output/Top_abstracts/stm17_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm17, n = 15)$prob)),path="Output/Top_words/stm17_words.csv")

#18 topics*******
number_of_topics=18
stm18thoughts0<-findThoughts(stm18,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm18thoughts<-findThoughts(stm18,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm18thoughts<-out$meta[data.frame(unlist(stm18thoughts))[,1],]
stm18thoughts$topic <-sort(rep(1:number_of_topics,25))

stm18theta <-data.frame(stm18$theta[data.frame(unlist(stm18thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm18theta)[1]){
  doc.thetas<-stm18theta[i,stm18thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm18thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm18_words<-data.frame((labelTopics(stm18, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm18_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm18thoughts=merge(stm18thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm18thoughts=stm18thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm18thoughts,path="Output/Top_abstracts/stm18_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm18, n = 15)$prob)),path="Output/Top_words/stm18_words.csv")

#19 topics*******
number_of_topics=19
stm19thoughts0<-findThoughts(stm19,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm19thoughts<-findThoughts(stm19,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm19thoughts<-out$meta[data.frame(unlist(stm19thoughts))[,1],]
stm19thoughts$topic <-sort(rep(1:number_of_topics,25))

stm19theta <-data.frame(stm19$theta[data.frame(unlist(stm19thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm19theta)[1]){
  doc.thetas<-stm19theta[i,stm19thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm19thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm19_words<-data.frame((labelTopics(stm19, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm19_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm19thoughts=merge(stm19thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm19thoughts=stm19thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm19thoughts,path="Output/Top_abstracts/stm19_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm19, n = 15)$prob)),path="Output/Top_words/stm19_words.csv")

#20 topics*******
number_of_topics=20
stm20thoughts0<-findThoughts(stm20,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm20thoughts<-findThoughts(stm20,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm20thoughts<-out$meta[data.frame(unlist(stm20thoughts))[,1],]
stm20thoughts$topic <-sort(rep(1:number_of_topics,25))

stm20theta <-data.frame(stm20$theta[data.frame(unlist(stm20thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm20theta)[1]){
  doc.thetas<-stm20theta[i,stm20thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm20thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm20_words<-data.frame((labelTopics(stm20, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm20_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))


stm20thoughts=merge(stm20thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm20thoughts=stm20thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm20thoughts,path="Output/Top_abstracts/stm20_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm20, n = 15)$prob)),path="Output/Top_words/stm20_words.csv")

#print topics with top words
stm20.top.words<-data.frame((labelTopics(stm20, n = 15)$prob))
stm20.top.words$theta.allyrs<-colMeans(stm20$theta[, 1:20]) #takes mean prevalence across the total corpus; theta is topic proportion matrix with document rows, topic columns
stm20.top.words$theta.recent<-colMeans(stm20$theta[which(out$meta$year>2015), 1:20]) #topic proportions, ordered by number 2016-19 YEARS ONLY

combine_words$theta.allyrs<-colMeans(stm20$theta[, 1:20])
combine_words$theta.recent<-colMeans(stm20$theta[which(out$meta$year>2015), 1:20]) #topic proportions, ordered by number 2016-19 YEARS ONLY
combine_words$theta.pre2016<-colMeans(stm20$theta[which(out$meta$year<=2015), 1:20]) #topic proportions, ordered by number 2016-19 YEARS ONLY

rev(order(colMeans(stm20$theta[, 1:20])))

#write_excel_csv(combine_words,path="Output/Top_words/stm20 topic categories.csv")

#21 topics*******
number_of_topics=21
stm21thoughts0<-findThoughts(stm21,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm21thoughts<-findThoughts(stm21,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm21thoughts<-out$meta[data.frame(unlist(stm21thoughts))[,1],]
stm21thoughts$topic <-sort(rep(1:number_of_topics,25))

stm21theta <-data.frame(stm21$theta[data.frame(unlist(stm21thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm21theta)[1]){
  doc.thetas<-stm21theta[i,stm21thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm21thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm21_words<-data.frame((labelTopics(stm21, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm21_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm21thoughts=merge(stm21thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm21thoughts=stm21thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm21thoughts,path="Output/Top_abstracts/stm21_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm21, n = 15)$prob)),path="Output/Top_words/stm21_words.csv")

#22 topics*******
number_of_topics=22
stm22thoughts0<-findThoughts(stm22,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm22thoughts<-findThoughts(stm22,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm22thoughts<-out$meta[data.frame(unlist(stm22thoughts))[,1],]
stm22thoughts$topic <-sort(rep(1:number_of_topics,25))

stm22theta <-data.frame(stm22$theta[data.frame(unlist(stm22thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm22theta)[1]){
  doc.thetas<-stm22theta[i,stm22thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm22thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm22_words<-data.frame((labelTopics(stm22, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm22_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm22thoughts=merge(stm22thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm22thoughts=stm22thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm22thoughts,path="Output/Top_abstracts/stm22_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm22, n = 15)$prob)),path="Output/Top_words/stm22_words.csv")

#23 topics*******
number_of_topics=23
stm23thoughts0<-findThoughts(stm23,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm23thoughts<-findThoughts(stm23,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm23thoughts<-out$meta[data.frame(unlist(stm23thoughts))[,1],]
stm23thoughts$topic <-sort(rep(1:number_of_topics,25))

stm23theta <-data.frame(stm23$theta[data.frame(unlist(stm23thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm23theta)[1]){
  doc.thetas<-stm23theta[i,stm23thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm23thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm23_words<-data.frame((labelTopics(stm23, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm23_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm23thoughts=merge(stm23thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm23thoughts=stm23thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm23thoughts,path="Output/Top_abstracts/stm23_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm23, n = 15)$prob)),path="Output/Top_words/stm23_words.csv")

#24 topics*******
number_of_topics=24
stm24thoughts0<-findThoughts(stm24,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm24thoughts<-findThoughts(stm24,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm24thoughts<-out$meta[data.frame(unlist(stm24thoughts))[,1],]
stm24thoughts$topic <-sort(rep(1:number_of_topics,25))

stm24theta <-data.frame(stm24$theta[data.frame(unlist(stm24thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm24theta)[1]){
  doc.thetas<-stm24theta[i,stm24thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm24thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm24_words<-data.frame((labelTopics(stm24, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm24_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm24thoughts=merge(stm24thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm24thoughts=stm24thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm24thoughts,path="Output/Top_abstracts/stm24_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm24, n = 15)$prob)),path="Output/Top_words/stm24_words.csv")

#25 topics (note that these model results are presented in main text of manuscript)
number_of_topics=25
stm25thoughts0<-findThoughts(stm25,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm25thoughts<-findThoughts(stm25,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm25thoughts<-out$meta[data.frame(unlist(stm25thoughts))[,1],]
stm25thoughts$topic <-sort(rep(1:number_of_topics,25))

stm25theta <-data.frame(stm25$theta[data.frame(unlist(stm25thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm25theta)[1]){
  doc.thetas<-stm25theta[i,stm25thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm25thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm25_words<-data.frame((labelTopics(stm25, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm25_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm25thoughts=merge(stm25thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm25thoughts=stm25thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm25thoughts,path="Output/Top_abstracts/stm25_top_abstracts.csv") #write to file

#print topics with top words
stm25.top.words<-data.frame((labelTopics(stm25, n = 15)$prob))
stm25.top.words$theta.allyrs<-colMeans(stm25$theta[, 1:25]) #takes mean prevalence across the total corpus; theta is topic proportion matrix with document rows, topic columns
stm25.top.words$theta.recent<-colMeans(stm25$theta[which(out$meta$year>2015), 1:25]) #topic proportions, ordered by number 2016-19 YEARS ONLY (for comparison)
#write_excel_csv(stm25.top.words,path="Output/Top_words/stm25_words.csv")

#26 topics*******
number_of_topics=26
stm26thoughts0<-findThoughts(stm26,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm26thoughts<-findThoughts(stm26,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm26thoughts<-out$meta[data.frame(unlist(stm26thoughts))[,1],]
stm26thoughts$topic <-sort(rep(1:number_of_topics,25))

stm26theta <-data.frame(stm26$theta[data.frame(unlist(stm26thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm26theta)[1]){
  doc.thetas<-stm26theta[i,stm26thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm26thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm26_words<-data.frame((labelTopics(stm26, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm26_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm26thoughts=merge(stm26thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm26thoughts=stm26thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm26thoughts,path="Output/Top_abstracts/stm26_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm26, n = 15)$prob)),path="Output/Top_words/stm26_words.csv")

#27 topics*******
number_of_topics=27
stm27thoughts0<-findThoughts(stm27,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm27thoughts<-findThoughts(stm27,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm27thoughts<-out$meta[data.frame(unlist(stm27thoughts))[,1],]
stm27thoughts$topic <-sort(rep(1:number_of_topics,25))

stm27theta <-data.frame(stm27$theta[data.frame(unlist(stm27thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm27theta)[1]){
  doc.thetas<-stm27theta[i,stm27thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm27thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm27_words<-data.frame((labelTopics(stm27, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm27_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm27thoughts=merge(stm27thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm27thoughts=stm27thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm27thoughts,path="Output/Top_abstracts/stm27_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm27, n = 15)$prob)),path="Output/Top_words/stm27_words.csv")

#28 topics*******
number_of_topics=28
stm28thoughts0<-findThoughts(stm28,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm28thoughts<-findThoughts(stm28,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm28thoughts<-out$meta[data.frame(unlist(stm28thoughts))[,1],]
stm28thoughts$topic <-sort(rep(1:number_of_topics,25))

stm28theta <-data.frame(stm28$theta[data.frame(unlist(stm28thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm28theta)[1]){
  doc.thetas<-stm28theta[i,stm28thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm28thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm28_words<-data.frame((labelTopics(stm28, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm28_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm28thoughts=merge(stm28thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm28thoughts=stm28thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm28thoughts,path="Output/Top_abstracts/stm28_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm28, n = 15)$prob)),path="Output/Top_words/stm28_words.csv")

#29 topics*******
number_of_topics=29
stm29thoughts0<-findThoughts(stm29,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm29thoughts<-findThoughts(stm29,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm29thoughts<-out$meta[data.frame(unlist(stm29thoughts))[,1],]
stm29thoughts$topic <-sort(rep(1:number_of_topics,25))

stm29theta <-data.frame(stm29$theta[data.frame(unlist(stm29thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm29theta)[1]){
  doc.thetas<-stm29theta[i,stm29thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm29thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm29_words<-data.frame((labelTopics(stm29, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm29_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm29thoughts=merge(stm29thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm29thoughts=stm29thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm29thoughts,path="Output/Top_abstracts/stm29_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm29, n = 15)$prob)),path="Output/Top_words/stm29_words.csv")


#30 topics*******
number_of_topics=30
stm30thoughts0<-findThoughts(stm30,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm30thoughts<-findThoughts(stm30,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm30thoughts<-out$meta[data.frame(unlist(stm30thoughts))[,1],]
stm30thoughts$topic <-sort(rep(1:number_of_topics,25))

stm30theta <-data.frame(stm30$theta[data.frame(unlist(stm30thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm30theta)[1]){
  doc.thetas<-stm30theta[i,stm30thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm30thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm30_words<-data.frame((labelTopics(stm30, n = 10)$prob)) #get top 10 words per topic
combine_words=unite(stm30_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm30thoughts=merge(stm30thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm30thoughts=stm30thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm30thoughts,path="Output/Top_abstracts/stm30_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm30, n = 15)$prob)),path="Output/Top_words/stm30_words.csv")

#31 topics*******
number_of_topics=31
stm31thoughts0<-findThoughts(stm31,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm31thoughts<-findThoughts(stm31,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm31thoughts<-out$meta[data.frame(unlist(stm31thoughts))[,1],]
stm31thoughts$topic <-sort(rep(1:number_of_topics,25))

stm31theta <-data.frame(stm31$theta[data.frame(unlist(stm31thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm31theta)[1]){
  doc.thetas<-stm31theta[i,stm31thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm31thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm31_words<-data.frame((labelTopics(stm31, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm31_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm31thoughts=merge(stm31thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm31thoughts=stm31thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm31thoughts,path="Output/Top_abstracts/stm31_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm31, n = 15)$prob)),path="Output/Top_words/stm31_words.csv")

#32 topics*******
number_of_topics=32
stm32thoughts0<-findThoughts(stm32,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm32thoughts<-findThoughts(stm32,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm32thoughts<-out$meta[data.frame(unlist(stm32thoughts))[,1],]
stm32thoughts$topic <-sort(rep(1:number_of_topics,25))

stm32theta <-data.frame(stm32$theta[data.frame(unlist(stm32thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm32theta)[1]){
  doc.thetas<-stm32theta[i,stm32thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm32thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm32_words<-data.frame((labelTopics(stm32, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm32_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm32thoughts=merge(stm32thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm32thoughts=stm32thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm32thoughts,path="Output/Top_abstracts/stm32_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm32, n = 15)$prob)),path="Output/Top_words/stm32_words.csv")

#33 topics*******
number_of_topics=33
stm33thoughts0<-findThoughts(stm33,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm33thoughts<-findThoughts(stm33,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm33thoughts<-out$meta[data.frame(unlist(stm33thoughts))[,1],]
stm33thoughts$topic <-sort(rep(1:number_of_topics,25))

stm33theta <-data.frame(stm33$theta[data.frame(unlist(stm33thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm33theta)[1]){
  doc.thetas<-stm33theta[i,stm33thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm33thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm33_words<-data.frame((labelTopics(stm33, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm33_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm33thoughts=merge(stm33thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm33thoughts=stm33thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm33thoughts,path="Output/Top_abstracts/stm33_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm33, n = 15)$prob)),path="Output/Top_words/stm33_words.csv")

#34 topics*******
number_of_topics=34
stm34thoughts0<-findThoughts(stm34,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm34thoughts<-findThoughts(stm34,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm34thoughts<-out$meta[data.frame(unlist(stm34thoughts))[,1],]
stm34thoughts$topic <-sort(rep(1:number_of_topics,25))

stm34theta <-data.frame(stm34$theta[data.frame(unlist(stm34thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm34theta)[1]){
  doc.thetas<-stm34theta[i,stm34thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm34thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm34_words<-data.frame((labelTopics(stm34, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm34_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm34thoughts=merge(stm34thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm34thoughts=stm34thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm34thoughts,path="Output/Top_abstracts/stm34_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm34, n = 15)$prob)),path="Output/Top_words/stm34_words.csv")

#35 topics*******
number_of_topics=35
stm35thoughts0<-findThoughts(stm35,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm35thoughts<-findThoughts(stm35,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm35thoughts<-out$meta[data.frame(unlist(stm35thoughts))[,1],]
stm35thoughts$topic <-sort(rep(1:number_of_topics,25))

stm35theta <-data.frame(stm35$theta[data.frame(unlist(stm35thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm35theta)[1]){
  doc.thetas<-stm35theta[i,stm35thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm35thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm35_words<-data.frame((labelTopics(stm35, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm35_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm35thoughts=merge(stm35thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm35thoughts=stm35thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm35thoughts,path="Output/Top_abstracts/stm35_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm35, n = 15)$prob)),path="Output/Top_words/stm35_words.csv")

#40 topics*******
number_of_topics=40
stm40thoughts0<-findThoughts(stm40,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm40thoughts<-findThoughts(stm40,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm40thoughts<-out$meta[data.frame(unlist(stm40thoughts))[,1],]
stm40thoughts$topic <-sort(rep(1:number_of_topics,25))

stm40theta <-data.frame(stm40$theta[data.frame(unlist(stm40thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm40theta)[1]){
  doc.thetas<-stm40theta[i,stm40thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm40thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm40_words<-data.frame((labelTopics(stm40, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm40_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm40thoughts=merge(stm40thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm40thoughts=stm40thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm40thoughts,path="Output/Top_abstracts/stm40_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm40, n = 15)$prob)),path="Output/Top_words/stm40_words.csv")

#50 topics*******
number_of_topics=50
stm50thoughts0<-findThoughts(stm50,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm50thoughts<-findThoughts(stm50,n = 25, texts=out$meta$Abstract.Title.Keywords,topics = c(1:number_of_topics))$index
stm50thoughts<-out$meta[data.frame(unlist(stm50thoughts))[,1],]
stm50thoughts$topic <-sort(rep(1:number_of_topics,25))

stm50theta <-data.frame(stm50$theta[data.frame(unlist(stm50thoughts0))[,1],])  #dataframe of relevant topic proportons from findThoughts (top 25 papers per topic),  rows are documents, columns are topics (now includes ALL topics, not just one topic per doc)

#for loop to pull out thetas for each top paper per topic
#define empty dataframes
doc.thetas <-NULL
doc.thetas2<-NULL
for(i in 1:dim(stm50theta)[1]){
  doc.thetas<-stm50theta[i,stm50thoughts$topic[i]];
  doc.thetas2<-rbind(doc.thetas2,doc.thetas)
}
stm50thoughts$theta<-round(doc.thetas2,digits=2) #add column of document level thetas (proportion of abstract assigned to topic), round to two decimals

stm50_words<-data.frame((labelTopics(stm50, n = 15)$prob)) #get top 15 words per topic
combine_words=unite(stm50_words,"top_words",sep=",")
combine_words$topic<-as.integer(row.names(combine_words))

stm50thoughts=merge(stm50thoughts,combine_words,by.x="topic",by.y="topic",all.x=TRUE) #append top words as column to top abstracts
stm50thoughts=stm50thoughts[c("id","authors","year","source","keywords","topics","Abstract.Title","theta","topic","top_words")] #remove some columns
write_excel_csv(stm50thoughts,path="Output/Top_abstracts/stm50_top_abstracts.csv") #write to file
#print topics with top words
#write_excel_csv(data.frame((labelTopics(stm50, n = 15)$prob)),path="Output/Top_words/stm50_words.csv")

#5) VISUALIZE model results (Fig. 3)-----------------------------

#read in custom functions for figures (Fig. 3)
source("Scripts/GBIF_literatureReview_customFunctions.R")

#plot topics by prevalence; expected topic proportion is the expected frequency across the entire corpus (the mean probability across all documents)
#mean proportions (=prevalence) across the total corpus; theta is topic proportion matrix with document rows, topic columns

#**25 topic model (presented in main text of paper):

#***note that topic numbers are NOT the same as number in paper (topics are numbered by descending frequency in paper for easier reader interpretation)
#note that "topic 12" is not all that meaningful -- removed from interpretation and figures (only accounts for 2.6% of corpus) (common thread seems to be plants, but little interpretive value)

# theta is topic proportion matrix with document rows, topic columns
colMeans(stm25$theta[, 1:25]) #topic proportions (entire corpus), vector is in order by number by model output;  
rev(order(colMeans(stm25$theta[, 1:25]))) #now in descending order of prevalence, number is model assigned topic number (e.g., first number is "topic 1" from model, number is rank order by topic prevalence)

frequency.1 <- colMeans(stm25$theta[, 1:25]) #topic proportions (entire corpus), vector is in order by number by model output; 
frequency <- colMeans(stm25$theta[, -c(12)]) #same as above, but remove "topic 12" for figures ("junk" topic)
frequency.focal<-colMeans(stm25$theta[which(out$meta$year>2015), 1:25]) #recent  2016-19 only

sort(frequency.focal,decreasing=T)
sum(frequency.focal[c(16,20)]) #SDM topics account for 13% in RECENT
sum(frequency.1[c(16,20)])  #SDM topics account for 13% in ALL YEARS

#vectors for topic names, ordered by model output (not ordered by prevalence as plotted in figures):
#Give meaningful topic names (all topics)
topicNames<-c("Population genetics","Biodiversity informatics","Novel occurrences/Range extensions","Invasion biology","Disease","Marine","Conservation","Regional distribution patterns",
              "Macroecoogical patterns", "Historical biogeography","Spatial ecology","Plant (***junk?)","Niche dynamics","Invasive species management","Ecophysiology","SDM (applications)",
              "Taxonomic descriptions and revisions","Global range dynamics","Clade diversification","SDM (tools,methods,theory)","Climate futures","People and nature (ethnobiology)","Forest biology","Natural history","Species interactions")

#Give meaningful topic names (all topics) excluding "junk" topic
topicNames.3<-c("Population genetics","Biodiversity informatics","Novel occurrences and range extensions","Invasion biology","Disease","Marine","Conservation","Regional distribution patterns",
                "Macroecological diversity patterns", "Historical biogeography","Spatial ecology","Niche dynamics","Invasive species management","Functional ecology","Species distribution models (applications)",
                "Taxonomic treatments","Global invasion dynamics","Clade diversification","Species distribution models (tools,theory)","Climate futures","People and nature","Forest biology","Phenotype","Species interactions")

#Barchart of topic proportions (no color codings)
par(mfrow=c(1,1),oma=c(0,0,0,0),xpd=F)
par(bty="n",col="black",lwd=5)
plot(stm25,type="summary",text.cex=0.9,n=10,xlim=c(0,.2),main=NA,xlab=NA,lwd=10,xaxt="n",topic.names=topicNames.3,topics=c(1:11,13:25),custom.labels = "")
axis(side=1, at=seq(0,0.15, by=0.05),lab=paste((seq(0,0.15, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)

#focal recent years only (2016-19); very similar to all years
plot.STM.focal(stm25,type="summary",text.cex=0.9,n=10,xlim=c(0,.2),main=NA,xlab=NA,lwd=10,xaxt="n",topic.names=topicNames.3,topics=c(1:11,13:25),custom.labels = "")
axis(side=1, at=seq(0,0.15, by=0.05),lab=paste((seq(0,0.15, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)

#Change color of lines to correspond to temporal change within topics (comparing recent 2016-19 studies to pre-2016 studies)

#calculate topic proportions for pre-2016 years  (to compare to recent slice of dataset, 2016-19)
frequency.pre2016<-colMeans(stm25$theta[which(out$meta$year<2016), 1:25])

#calculate topic change as a proportion within topics
topicChanges<-(frequency.focal-frequency.pre2016)/frequency.pre2016
topicChanges<-topicChanges[-c(12)] #remove "junk" topic

#load packages for color gradient
library(RColorBrewer)
library(scales)
library(plotrix) #for legend

colors <- brewer.pal(8, "RdYlBu") #RdYlBu
pal <- colorRampPalette(colors)
heat.palette = rev(pal(1001)) #1001 equally spaced colors
#show_col(heat.palette[500]) #light tan is no change

heat.palette.neg = heat.palette[499:1] #subset vector for negative colors (blues)
heat.palette.pos = heat.palette[500:1001] #subset vector for positive colors (reds)

topicColors = rep(NA,length(topicChanges)) #empty vector for colors showing % change between pre and post 2016
topicColors[topicChanges<0] = heat.palette.neg[round((topicChanges[topicChanges<0]*-1)/max(topicChanges)*500)] #more blue is more negative
topicColors[topicChanges>0] = heat.palette.pos[round((topicChanges[topicChanges>0])/max(topicChanges)*500)] #more red is more positive

topicColors = topicColors[order(frequency)] #reorder by increasing topic prevalence (overall)

par(mfrow=c(1,1),oma=c(0,0,0,0),xpd=F)
par(bty="n",col="black",lwd=10)
plot.STM.width(stm25,type="summary",text.cex=0.9,n=10,xlim=c(0,.2),main=NA,xlab=NA,lwd=10,xaxt="n",topic.names=topicNames.3,topics=c(1:11,13:25),custom.labels = "")
axis(side=1, at=seq(0,0.15, by=0.05),lab=paste((seq(0,0.15, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)
par(bty="n",col="black",lwd=1)
color.legend(0.14,0.9,0.15,6.3,rect.col=heat.palette,gradient="y",align="rb",cex=1)#legend=seq(-0.47,0.47,length=5)

#version with no labels (to add different labels or improve aesthetics in editing programs outside of R)
par(mfrow=c(1,1),oma=c(0,0,0,0),xpd=F)
par(bty="n",col="black",lwd=10)
plot.STM.width(stm25,type="summary",text.cex=0.9,n=10,xlim=c(0,.2),main=NA,xlab=NA,lwd=10,xaxt="n",topics=c(1:11,13:25),custom.labels = "",topic.names="")
axis(side=1, at=seq(0,0.15, by=0.05),lab=paste((seq(0,0.15, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)
par(bty="n",col="black",lwd=1)
color.legend(0.14,0.9,0.15,6.3,rect.col=heat.palette,gradient="y",align="rb",cex=1)

#topic correlation graph
stm25.network <- topicCorr(stm25)

plot(stm25.network,vertex.color="white",vertex.shape="fcircle",vertex.frame.width=5,vertex.size=500*frequency,edge.width=2,edge.color="darkgrey",vertex.label.dist=0,vertex.label.cex=0.75,vertex.label.family="sans",edge.lty=1,vlabels=topicNames.3,vertex.label.font=2,topics=c(1:11,13:25),vertex.frame.color=topicColors)
#version with no labels for editing outside R:
plot(stm25.network,vertex.color="white",vertex.shape="fcircle",vertex.frame.width=5,vertex.size=500*frequency,edge.width=2,edge.color="darkgrey",vertex.label.dist=0,vertex.label.cex=0.75,vertex.label.family="sans",edge.lty=1,vlabels=NA,vertex.label.font=2,topics=c(1:11,13:25),vertex.frame.color=topicColors)
#add color legend
color.legend(-1.2,-0.73,-1.5,-1.3,rect.col=heat.palette,gradient="y",align="rb",cex=1)#legend=seq(-0.47,0.47,length=5)

#Plots stm 20 (supplemental figure for comparison)--------------

#Give meaningful topic names (all topics), "junk" or different from 25 topic model denoted with "*"
topicNames.20<-c("Forest biology","Biodiversity informatics","Taxonomy+Novel occurrences*","Invasion biology","Disease+Species interactions*","Climate change+Marine*","Conservation","Species distributions","Macroecological patterns","Populations genetics","Land use+SDMs*","Ethnobotany*","SDMs (applications)","Invasive species management","Functional ecology","Climate futures","Taxonomic treatments+Phylogeography*","Global invasion dynamics","Historical biogeography","Species distributions II+Range extensions*")

#Barchart of topic proportions
par(mfrow=c(1,1),oma=c(0,0,0,0),xpd=F)
par(bty="n",col="black",lwd=5)

plot(stm20,type="summary",text.cex=0.9,n=10,xlim=c(0,.25),main=NA,xlab=NA,lwd=10,xaxt="n")
plot(stm20,type="summary",text.cex=0.9,n=10,xlim=c(0,.20),main=NA,xlab=NA,lwd=10,xaxt="n",topic.names=topicNames.20,custom.labels = "")
axis(side=1, at=seq(0,0.20, by=0.05),lab=paste((seq(0,0.20, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)

#Plots stm 30 (supplemental figure for comparison)-----------------

#Give meaningful topic names (all topics), "junk" or different from 25 topic model denoted with "*"
topicNames.30<-c("Forest biology","Biodiversity informatics","Range extensions*","Invasion biology","Crop pests*","Marine","Conservation","Distribution patterns I*","Macroecological patterns","Historical biogeography","Distribution patterns II+Conservation*","Botany (?)*","Niche dynamics","Invasive species management","Functional ecology","SDMs (applications)","Taxonomy I*","Global invasion dynamics","Taxonomy II*","SDMs (methods,tools,theory)","Climate futures","Ethnobotany*","Habitat+Conservation+Landscape ecology*","Phenotype+Bird phenology*","Evolution (polyploidy,hybridization,gene flow)","Clade diversification","Population genetics","Novel occurrences*","Disease","Thermal biology*")

#Barchart of topic proportions
par(mfrow=c(1,1),oma=c(0,0,0,0),xpd=F)
par(bty="n",col="black",lwd=5)

plot(stm30,type="summary",text.cex=0.9,n=10,xlim=c(0,.25),main=NA,xlab=NA,lwd=10,xaxt="n")
plot(stm30,type="summary",text.cex=0.9,n=10,xlim=c(0,.30),main=NA,xlab=NA,lwd=10,xaxt="n",topic.names=topicNames.30,custom.labels = "")
axis(side=1, at=seq(0,0.15, by=0.05),lab=paste((seq(0,0.15, by=0.05))*100,"%"),cex.axis=1.5)
mtext(expression("Topic Proportions"),side=1,cex=1.8,line=3.5)

