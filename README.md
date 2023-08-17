# MGRec
**Cristian Díaz-Muñoz** and **Marko Verce**

*The actual plots that result from running the codes below can't be shown, since this dataset is still unpublished. The genus and species name of the case study is shown as "Genus species" to mantain the confidentiality.*

## Introduction

This is a guide to perform metagenomic reads recruitment plotting. Most of the code has been developed by Marko Verce, with some modifications by me. You will need two files to start this script:

  1. **MGrec.xxx.besthit.tab**. A tab separated file containing the best hit for each read against a BLAST database.
  2. **xxx_contig_len_gen_sp.txt**. A tab separated file containing the length, genus name, and species name for every contig. 

The information on how to generate those files is publicly available in Verce *et al.,* 2019 (10.3389/fmicb.2019.00479).

## Usage

This script is based on `MGRec3_graphs.R` adapted for a reduced number of samples. Most of the code is original of Marko Verce, so if you want to acknowledge this work, please do so mostly to him. In this guide we will use samples from a colleague, who was analyzing the microbiota of some spontaneous cocoa fermentation processes in a tropical country. In this case, a single metagenomic read file was aligned toward a database containing the genomes of type strains for every species contained in 25 different genera of bacteria and fungi. 

First of all, let's set the working directory and the load the libraries needed for this script.

```{r}
setwd("C:/path/to/your/working/directory")

library(scales) # Scale functions for visualization
library(tidyverse) # Type of R syntax, to make code easier
library(grid) # Add grid to a plot
library(gridExtra) # Miscellaneous functions for "grid" graphics
library(ggpubr) # To plot multiple objects together
```


### Parsing

We are going to load functions to the workspace to read and obtain the information that we want in the format that we want from the two input files. Let's start with a function that will be used later by another function. It basically takes the length of all contigs within a species (`lengths`), calculate the cumulative sum of all contigs (`cumlengths`) and then calculates the starting point of each contig (`cumsum0`) so that each one gets ordinated (plotted) right after the next one. 

```{r}
cumsum0<-function(lengths){
  cumlengths<-cumsum(as.numeric(lengths))
  cumsum0<-cumlengths-lengths
  return(cumsum0)
}
```

The following function takes the list of contigs (accession, length, genus, species) used in the blast database and adds numbers used to help create the X coordinates plotting the species. Also, we use the function `mutate()` to order the tibble. This is very necessary, because `tapply` is not intrinsically connected to the tibble and apparently R follows a different alphabetic order than Python (Which is first? *Bacillus soli* or *Bacillus solimangrovi*?).Due to the ordering difference, the `Offset_sp` vector is crap unless you order the tibble to R order.

```{r}
List_sp_Offset<-function(file){
  A <- read_tsv(file,col_names = c("Accession","Length","Genus","Species")) %>%
    arrange(Genus,Species) %>% 
    mutate(Offset_sp = unlist(tapply(.$Length,.$Species,cumsum0),use.names=F))
  return(A)
}
```

The following function takes the output of `List_sp_Offset` and is used to calculate offsets for genus plots.

```{r}
List_ge_Offset<-function(df){
  B <- df %>%
    arrange(Genus,Species) %>%
    group_by(Species,Genus) %>% 
    summarise(Length=sum(Length)) %>% 
    ungroup() %>%
    mutate(Offset_ge_basis=unlist(tapply(.$Length,.$Genus,cumsum0),use.names=F),
           Cumsum=unlist(tapply(.$Length,.$Genus,function(x) cumsum(as.numeric(x))),use.names=F))
  return(B)
}
```

The following function reads the filtered blast output, uses the contig/species lists with genus/species offsets, and creates the X coordinates for genus-level and species-level figures.

```{r}
processBLASTbest_memsave<-function(file,List_sp_Offset,List_ge_Offset,qcov){
  BL<-read_tsv(file,col_names=c("query","Accession","pid","len","mismatch",
                                "gapopen","qstart","qend","sstart","send",
                                "evalue","bitscore","qlen","slen")) %>%
    left_join(List_sp_Offset,by="Accession") %>%
    left_join(List_ge_Offset,by="Species") %>%
    mutate(X_sp = (sstart+send)/2 + Offset_sp,
           X_ge = (sstart+send)/2 + Offset_ge_basis + Offset_sp,
           qcovmod=(qend-qstart+1)/qlen*100) %>%
    select(Query=query,Accession=Accession,Percentage_id=pid,Length=len,
           Mismatch=mismatch,GapOpen=gapopen,Query_start=qstart,Query_end=qend,
           Subject_start=sstart,Subject_end=send,Evalue=evalue,Bitscore=bitscore,
           Query_length=qlen,Subject_length=slen,Query_coverage=qcovmod,
           Genus=Genus.x,Species=Species,Offset_species=Offset_sp,
           Offset_genus_basis=Offset_ge_basis,X_species=X_sp,X_genus=X_ge) %>%
    select(-Length,-Mismatch,-GapOpen,-Query_start,-Query_end,-Offset_species,
           -Offset_genus_basis,-Evalue,-Bitscore) %>%
    filter(Query_coverage>qcov)
  return(BL)
}
```

List of contig names, lengths (of the contigs), genus, species, and the cumulative sum of the contigs per species:

```{r}
Acc_len_gen_sp <- List_sp_Offset("XXX_contig_len_gen_sp.txt")
head(Acc_len_gen_sp)
# If we want to save this file as a text file to inspect it on Excel, for example, 
# we can save it. Also, it will be useful to restart from here next time.
write_tsv(Acc_len_gen_sp,"Acc_len_gen_sp_XXX.txt",col_names=F)
```

List the species, genus, lengths (of the species genomes), the cumulative sum (0) and the cumulative sum:

```{r}
Sp_gen_Genome <- List_ge_Offset(Acc_len_gen_sp)
head(Sp_gen_Genome)
# If we want to save this file as a text file to inspect it on Excel, for example, 
# we can save it. Also, it will be useful to restart from here next time.
write_tsv(Sp_gen_Genome,"Sp_gen_Genome_XXX.txt",col_names=F)
```

Import the BLAST results. As we saw above, this function takes the tabular file output of BLAST with only one hit per query sequence (read), the two dataframes that we just imported with the lengths of the contigs and species, and the minimum query cover to consider an aligment (any alignment below that number will be filtered out). Typically, this minimum query cover number is the same as the one we set during the BLAST alignment.

```{r}
df <- processBLASTbest_memsave("MGrec.XXX.XXX.besthit.tab",Acc_len_gen_sp, Sp_gen_Genome,60)
head(df)
```

### Generate tables

Now, we are going to generate the table for each species or genus and their percentage of reads recruited. The only valid interpretation is: "This species/genus recruited xx % of all reads that mapped against the database and survived the filters". This is because the database we constructed does not contain all microorganisms in the world (> 99 % of all microorganisms in the planet are not known), and depending on your matrix, there may be reads belonging to non-microbial DNA (as is the case for *Theobroma cacao* DNA in metagenomic samples of cocoa fermentation processes). So, to do that, we are using a function that calculates the number of unique reads that aligned to each species/genus (`taxlevel`) and divide it by the total amount of reads in the sample (`readnumber`). Finally, it arranges the dataframe by the percentage of reads recruited in descending order.

```{r}
sp_gen_table_filtered2<-function(df,taxlevel,readnumber){
  taxlevel <- enquo(taxlevel)
  df <- df %>% 
    group_by_at(vars(!!taxlevel)) %>% 
    summarise(Percentage=length(unique(Query))/readnumber*100) %>%
    ungroup() %>%
    arrange(desc(Percentage))
  return(df)
}
```

We can now run the function on our sample, in this case, using the species as the `taxlevel`:

```{r}
df_perc <- sp_gen_table_filtered2(df,"Species",nrow(df))
head(df_perc)
# If we want to save this file as a text file to inspect it on Excel, for example, 
# we can save it.
write.table(x = df_perc, file = "XXX_sp_recruitment_percentage.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE)
```


### Fragment recruitment plotting

The following function take the processed BLAST tibble, a species name, and makes a species-level plot. With `vircol` you can specify a continuous color palette from the `viridis()` package. The continuous scale with color the normalised counts of reads per locus of the genome. This is useful to see if there are some over represented loci. We define the `binwidth` in relation to the genome size of each species, as well as the `limits` of the x-axis.

```{r}
RecruitmentPlot_Species <- function(df, species, Acc_len_gen_sp, Sp_gen_Genome, vircol,
                                    binwidth.x=c(Sp_gen_Genome$Length[Sp_gen_Genome$Species==species]/1000),
                                    binwidth.y=0.25,
                                    limits=c(0,Sp_gen_Genome$Length[Sp_gen_Genome$Species==species]*1.001),
                                    full_y_axis=F,guide=T) {
  df <- filter(df,Species==species)
  Acc_len_gen_sp <- Acc_len_gen_sp[Acc_len_gen_sp$Species==species,]
  theme_base <- theme(axis.ticks.x=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x= element_blank(),
                      axis.line.y=element_line(colour="black"),
                      axis.ticks.y= if(full_y_axis==T) element_line(colour="black") else element_blank(),
                      axis.title.y= if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      axis.text.y= if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      panel.background=element_blank(),
                      panel.grid.major.x=element_blank(),
                      panel.grid.minor.x=element_line(colour="gray50"),
                      panel.grid.major.y=element_blank(),
                      panel.grid.minor.y=element_blank(),
                      plot.title=element_text(hjust=0.5))
  PLOT<- ggplot(df,aes(x=X_species,y=Percentage_id))+
    geom_bin2d(binwidth=c(binwidth.x,binwidth.y), aes(fill=..count../max(..count..)))+
    ggtitle(substitute(italic(b), list(b=species))) +
    theme_base+
    scale_y_continuous(labels=comma, limits=c(50,100.5), expand=c(0,0))+
    scale_x_continuous(name=NULL,breaks=NULL,
                       labels=NULL,
                       minor_breaks=NULL,
                       limits=limits,
                       expand=c(0,0))
  if (guide==T) {
    PLOT <- PLOT +
      scale_fill_viridis_c(name = "Normalised counts", 
                           option = vircol, end = 0.75, direction = -1)
  }
  else {
    PLOT <- PLOT +
      scale_fill_viridis_c(name = "Normalised counts", 
                           option = vircol, direction = -1, guide = F)
  }
  
  if (full_y_axis==T) {
    PLOT <- PLOT +
      scale_y_continuous("Identity [%]", labels=comma, 
                         limits=c(50,100.5), expand=c(0,0))
  }
  return(PLOT)
}
```

We can use this function to plot every species we are interested in. As we can be interested in many species, is better to save them as a variable, especially if we want to plot them together afterwards. In this example, I suspect that there is a putative new species close to *Genus species*: 

```{r}
a <- RecruitmentPlot_Species(df,"Genus species",Acc_len_gen_sp,Sp_gen_Genome, 
                             "cividis", full_y_axis = T,guide=T)
plot(a)
```

We see that the reads cover more or less the whole genome sequence of the species, albeit with a low coverage. We can deduce that it is probably low in relative abundance. If we look at the `df_perc` table we can see that *G. species* has recruited 0.73 % of all reads mapped. Further, we see that on average, the reads show an identity percentage around 90 %. Usually, strains of the same species share an identity equal or superior to 95 % and different species of the same genus show an identity around 80 %, so this is a quite controversial case. If you want to plot multiple species, from the same genus or different simoultaneously, we can save every plot in a different variable and plot them together using `garange()`:

```{r eval=FALSE}
# This is just an example where I want to plot 10 different species. 
# We can also control other parameters, as the legend, number of columns and 
# number of rows. If you want to control more aesthetics, type ?ggarange().

ggarrange(a,b,c,d,e,f,g,h,i,j,common.legend = TRUE,legend="right", ncol = 2, nrow = 5)
```


But now, let's say that we want to see all the *Genus* present in the sample. We could run this function for every species in the database, which would take us some time. Alternatively, we can also run a similar function that plots all species in a genus in a single plot. To do that, we use the `facet_wrap()` function of `ggplot`. If you set the number of rows for the facet as 0 ( by default it is `rowsize=0`), it will automatically decide the optimal number of rows for the plot, but if you're not happy with this decision, you can always set the number of rows yourself (*e.g..,* `rowsize=2`).

```{r}
RecruitmentPlot_Genus_CDM <- function(df, genus, Acc_len_gen_sp, Sp_gen_Genome, vircol,
                                      full_y_axis=F,guide=T,rowsize=0) {
  df <- filter(df,Genus==genus)
  g <- filter(Sp_gen_Genome, Genus==genus)
  h <- Acc_len_gen_sp[Acc_len_gen_sp$Genus==genus,]
  binwidth.x=c(mean(g$Length)/1000)
  binwidth.y=0.25
  if (rowsize==0) {
    if (nrow(g) <= 6) {
      rowsize <- 1
    }
    else if (nrow(g) > 6 & nrow(g) <= 12) {
      rowsize <- 2
    }
    else if (nrow(g) > 12 & nrow(g) <= 19) {
      rowsize <- 3
    }
    else if (nrow(g) >= 20 & nrow(g) < 30) {
      rowsize <- 4
    }
    else if (nrow(g) >= 30) {
      rowsize <- 5
    }
  }
  theme_base <- theme(axis.title=element_blank(),
                      axis.text.x= element_blank(),
                      axis.line.y=element_line(colour="black"),
                      axis.ticks.y= if(full_y_axis==T) element_line(colour="black") else element_blank(),
                      axis.title.y= if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      axis.text.y= if(full_y_axis==T) element_text(colour="black") else element_blank(),
                      panel.background=element_blank(),
                      panel.grid.major.x=element_blank(),
                      panel.grid.minor.x=element_line(colour="gray50"),
                      panel.grid.major.y=element_blank(),
                      panel.grid.minor.y=element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "italic"))
  PLOT<-ggplot(df,aes(x=X_species,y=Percentage_id))+
    geom_bin2d(binwidth=c(binwidth.x,binwidth.y), aes(fill=..count../max(..count..)))+
    geom_vline(data=g, aes(xintercept = Length)) +
    geom_vline(xintercept = 0) +
    facet_wrap(.~Species, nrow = rowsize, scales = "free") + 
    theme_base+
    scale_y_continuous(labels=comma, limits=c(60,100.5), expand=c(0,0))+
    scale_x_continuous(name=NULL,breaks=NULL,
                       labels=NULL,
                       minor_breaks=NULL,
                       expand=c(0,0))
  if (guide==T) {
    PLOT <- PLOT +
      scale_fill_viridis_c(name = "Normalised counts", option = vircol, end = 0.75, direction = -1)
  }
  else {
    PLOT <- PLOT +
      scale_fill_viridis_c(name = "Normalised counts", option = vircol, direction = -1, guide = F)
  }
  
  if (full_y_axis==T) {
    PLOT <- PLOT +
      scale_y_continuous("Identity [%]", labels=comma, limits=c(60,100.5), expand=c(0,0))
  }
  #return(rowsize)
  return(PLOT)
}
```

Okay, let's specify the *Genus* genus and, to change colors, we will use the "magma" palette this time. Let's set the number of rows to 2. Don't worry about the names, we can control the size of every element when saving the figure with `ggsave()`.

```{r}
b <- RecruitmentPlot_Genus_CDM(df,"Genus",Acc_len_gen_sp,Sp_gen_Genome, 
                               "magma", full_y_axis = T,guide=T, rowsize = 2)
plot(b)
```

We can see now that in the sample there is also *Genus species*, of whom we are certain about as most of the aligned reads are between 95 and 100 % identity. We also see that the coverage is much higher than for *G. species*, which correlates with the higher percentage of recruited reads (5.23 %). 

Finally, we are going to save the last plot. I typically save the plots in a different folder, in `tiff` format and with 300 dpi for optimal resolution suitable for publications. With `width` and `height` you can control the dimensions of the figure. This is pure trial and error until you get the figure as you want. A rule of thumb is that the higher the `width` and `height` are, the smaller will be the text elements in the text, but also the bigger will be the size (MB) of the file. 

```{r eval=FALSE}
ggsave("Plots/Genus.tiff", width=10,height=8, dpi=300)
```
