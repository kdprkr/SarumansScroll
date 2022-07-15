### ************************************
### PREFACE ----
### ************************************

# this script houses the entirety of code for analysis performed in R^
# ^after sequence processing performed with QIIME 2 (see Script01.R)
# ** note for KDP: PREFACE version 0.7 ** #

# configure the script
# NOTE: my machine runs MacOS; paths have not been tested on other platforms
# NOTE: this information and more is stored in section I below

setwd("~/Desktop/") # *gasp*
date_ofrun <- "(January 06, 2021 (01_06_21))" # date script was run
name_scrpt <- "Script02" # script filename
name_projd <- "SarumansScroll" # name of project directory
path_projd <- paste(sep = "", getwd(), "/", name_projd) # path to project dir
path_vault <- paste(sep = "", path_projd, "/vault") # path to storage vault
path_mtabs <- paste(sep = "", path_projd, "/mtabs") # path to metabolite dir
wksp_extsn <- paste(sep = "", name_scrpt, ".Rdata") # recyclable ext for outputs
text_extsn <- paste(sep = "", name_scrpt, ".txt") # recyclable ext for outputs
pdfs_extsn <- paste(sep = "", name_scrpt, ".pdf") # recyclable ext for outputs
PREFACE <- c("date_ofrun", "name_scrpt", "name_projd", "path_projd", 
             "path_vault", "path_mtabs", 
             "wksp_extsn", "text_extsn", "pdfs_extsn")

## INPUTS
# paths to universal inputs stored in the project directory (path_projd)
ifp_smp_dat <- paste(sep = "", path_projd, "/sampledata_SS.txt")

# paths to section specific inputs stored in the central vault (path_vault)
A.ifv_fts_tab <- paste(sep = "", path_vault, "/q2_fts_tab.tsv")
A.ifv_rep_seq <- paste(sep = "", path_vault, "/q2_rep_seq.fasta")
A.ifv_int_ggs <- paste(sep = "", path_vault, "/q2_int_ggs.tsv")
A.ifv_int_slv <- paste(sep = "", path_vault, "/q2_int_slv.tsv")
E.ifv_chao <- paste(sep = "", path_vault, "/q2_chao.tsv")
E.ifv_shan <- paste(sep = "", path_vault, "/q2_shan.tsv")

# paths to section specific inputs stored in the metabolite dir (path_mtabs)
Y.ifm_smp_dat_mtb <- paste(sep = "", path_mtabs, "/sampledata_mtabs_SS.txt")
Y.ifm_cvcln_mtb <- paste(sep = "", path_mtabs, "/data_mtabs_cvColon.txt")
Y.ifm_cvbld_mtb <- paste(sep = "", path_mtabs, "/data_mtabs_cvPlasma.txt")
Y.ifm_gfcln_mtb <- paste(sep = "", path_mtabs, "/data_mtabs_gfColon.txt")
Y.ifm_gfbld_mtb <- paste(sep = "", path_mtabs, "/data_mtabs_gfPlasma.txt")

## OUTPUTS
# paths for outputs stored in the central vault (path_vault)
ofv_COSMOS_wksp <- paste(sep = "", path_vault, "/WS_COSMOS_", wksp_extsn)

I.ofv_info <- paste(sep = "", path_vault, "/info_Section_I_", text_extsn)
I.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_I_", wksp_extsn)

A.ofv_info <- paste(sep = "", path_vault, "/info_Section_A_", text_extsn)
A.ofv_prov <- paste(sep = "", path_vault, "/prov_Section_A_", text_extsn)
A.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_A_", wksp_extsn)
A.ofv_EBTKS_raw <- paste(sep = "", path_vault, "/table_EBTKS_abs_raw.txt")
A.ofv_EBTKS_pro <- paste(sep = "", path_vault, "/table_EBTKS_abs_processed.txt")

S.ofv_info <- paste(sep = "", path_vault, "/info_Section_S_", text_extsn)
S.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_S_", wksp_extsn)

T.ofv_prov <- paste(sep = "", path_vault, "/prov_Section_T_", text_extsn)
T.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_T_", wksp_extsn)
T.ofv_EBTKS_pro_tree <- paste(sep = "", path_vault, "/tree_EBTKS_processed.nwk")

E.ofv_info <- paste(sep = "", path_vault, "/info_Section_E_", text_extsn)
E.ofv_prov <- paste(sep = "", path_vault, "/prov_Section_E_", text_extsn)
E.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_E_", wksp_extsn)
E.ofv_smp_adiv <- paste(sep = "", path_vault, "/table_adiv_with_sampledata.txt")
E.ofv_rst_adiv_EPD <- paste(sep = "", path_vault, "/table_adiv_all_stats.txt")
E.ofv_rst_flt_adiv <- paste(sep = "", path_vault, "/table_adiv_flt_stats.txt")
E.ofv_gpr_chao_cec <- paste(sep = "", path_vault, "/plot_adiv_chao_cecum.pdf")
E.ofv_gpr_shan_cec <- paste(sep = "", path_vault, "/plot_adiv_shan_cecum.pdf")
E.ofv_gpr_chao_cln <- paste(sep = "", path_vault, "/plot_adiv_chao_colon.pdf")
E.ofv_gpr_shan_cln <- paste(sep = "", path_vault, "/plot_adiv_shan_colon.pdf")
E.ofv_gga_cec <- paste(sep = "", path_vault, "/plot_adiv_cecum.pdf")
E.ofv_gga_cln <- paste(sep = "", path_vault, "/plot_adiv_colon.pdf")

B.ofv_prov <- paste(sep = "", path_vault, "/prov_Section_B_", text_extsn)
B.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_B_", wksp_extsn)
B.ofv_rst_bdiv_EPD <- paste(sep = "", path_vault, "/table_bdiv_all_stats.txt")
B.ofv_rst_flt_bdiv <- paste(sep = "", path_vault, "/table_bdiv_flt_stats.txt")
B.ofv_gpr_guni_cec <- paste(sep = "", path_vault, "/plot_bdiv_guni_cecum.pdf")
B.ofv_gpr_atch_cec <- paste(sep = "", path_vault, "/plot_bdiv_atch_cecum.pdf")
B.ofv_gpr_guni_cln <- paste(sep = "", path_vault, "/plot_bdiv_guni_colon.pdf")
B.ofv_gpr_atch_cln <- paste(sep = "", path_vault, "/plot_bdiv_atch_colon.pdf")
B.ofv_gga_cec <- paste(sep = "", path_vault, "/plot_bdiv_cecum.pdf")
B.ofv_gga_cln <- paste(sep = "", path_vault, "/plot_bdiv_colon.pdf")

Y.ofv_prov <- paste(sep = "", path_vault, "/prov_Section_Y_", text_extsn)
Y.ofv_wksp <- paste(sep = "", path_vault, "/WS_Section_Y_", wksp_extsn)
Y.ofv_gpr_cvcln_mtb <- paste(sep = "", path_vault, "/plot_pca_mtab_cvColon.pdf")
Y.ofv_gpr_cvbld_mtb <- paste(sep = "", path_vault, "/plot_pca_mtab_cvBlood.pdf")
Y.ofv_gpr_gfcln_mtb <- paste(sep = "", path_vault, "/plot_pca_mtab_gfColon.pdf")
Y.ofv_gpr_gfbld_mtb <- paste(sep = "", path_vault, "/plot_pca_mtab_gfBlood.pdf")

# save PREFACE workspace
PREFACE.lst <- c(ls(pattern = "ifp"), ls(pattern = "ifv"),ls(pattern = "ifm"),
                 ls(pattern = "ofv"), PREFACE, "PREFACE", "PREFACE.lst")
save(list = PREFACE.lst, 
     file = paste(sep = "", path_vault, "/WS_PREFACE_", wksp_extsn))

### ************************************
### UNIVERSAL OBJECTS ~ COSMOS ----
### ************************************

# this section contains items that are universal across datasets

# read in the sample data (commonly called metadata)
smp_dat <- read.table(ifp_smp_dat, sep = "\t", header = T, as.is = T, 
                      stringsAsFactors = F)

# hex codes for a gradient of black/grey/white (left to right = dark to light)
greydient <- c("#000000", "#252525", "#525252",
               "#969696", "#bbbbbb", "#d9d9d9", "#e0e0e0", 
               "#ffffff")

# purple gradient
purples <- c("#542788", "#8073ac", "#f2ebf9")

# Life Aquatic with Steve Zissou custom color palettes
ziss_blue <- c("#218ec4", "#2ca1db", "#67bbe5", "#92ceec", "#e9f5fb", "#ffffff")
ziss_navy <- c("#0f1b3e", "#192e67", "#2d52b9", "#466cd2", "#eaeffa", "#ffffff")
ziss_reds <- c("#b31300", "#f21a00", "#ff604d", "#ffbbb3", "#ffe8e5", "#ffffff")
ziss_ylow <- c("#fec72A", "#fede80", "#ffebb3", "#fff8e6", "#ffffff")
ziss_teal <- c("#285854", "#489e97", "#84c7c2", "#cae7e5", "#edf7f6", "#ffffff") 
ziss_grey <- c("#252525", "#3d3d3d", "#c4c4c4", "#d9d9d9", "#e9e9e9", "#ffffff")

# font family
# fnt <- "Courier"
fnt <- "Helvetica"

# multi-panel label ~ font parameters
pan_fnt <- list(size = 13, family = fnt, face = "bold")

# plot widths for single and double column
wid_sgl <- 78
wid_dbl <- 165

# vector with column names commonly removed or retained
com_col <- c("FeatureID", "RepSeq", 
             "int.ggs.tax", "int.ggs.Kingdom", "int.ggs.Phylum", 
             "int.ggs.Class", "int.ggs.Order", "int.ggs.Family", 
             "int.ggs.Genus", "int.ggs.Species", "int.ggs.cnf", 
             "int.ggs.lws.txn", "int.ggs.lws.lvl", 
             "int.slv.tax", "int.slv.Kingdom", "int.slv.Phylum", 
             "int.slv.Class", "int.slv.Order", "int.slv.Family", 
             "int.slv.Genus", "int.slv.Species", "int.slv.cnf", 
             "int.slv.lws.txn", "int.slv.lws.lvl")

# ^^cutpoints and symbols to be used for determining P & BH_P value significance
# NOTE: used as argumants for symnum() base function
cutpts <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
symbls <- c("****", "***", "**", "*", " ")

cutpts_P <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
symbls_P <- c("****", "***", "**", "*", " ")
cutpts_BHP <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
symbls_BHP <- c("****", "***", "**", "*", "ns")

# create a vector naming all of the above (and include ifp_smp_dat from PREFACE)
# NOTE: this is useful when saving workspaces
cosmos <- c("ziss_blue", "ziss_navy", "ziss_reds", "ziss_ylow", "ziss_teal", 
            "ziss_grey", "purples",
            "ifp_smp_dat", "smp_dat", "greydient", "fnt", "pan_fnt", "com_col",
            "cutpts", "symbls",
            "cutpts_P", "symbls_P", "cutpts_BHP", "symbls_BHP",
            "wid_sgl", "wid_dbl")

### ************************************
### ^^^UNIVERSAL OBJECTS ~ dataset/project specific ----
### ************************************

# this section contains items that are not universal across datasets ...
# ... but are universal across sections within a specific dataset or project

# ^^^abbreviations used below:
# _shp = shape; _lne = line
# _cec = cEcum; _prx = Proximal colon; _dtl = Distal colon; _fec = feces
# _C = Control; # _R = Rice_bran; _F = Fermented_rice_bran

# ggplot2 codes for shapes by Type and DietGroup
shp_cec <- 24 # open triangle (outline dark; fill dark or medium)
shp_prx <- 21 # open circle (outline medium; fill lightest)
shp_dtl <- 21 # open circle (outline light; fill light)
shp_fec <- 04 # x
# shp_C <- 67 # capital C
# shp_R <- 82 # capital R
# shp_F <- 70 # fapital F
shp_C <- 21 # fillable circle
shp_R <- 24 # fillable triangle (point at top)
shp_F <- 25 # fillable triangle (point at bottom)

# hex codes to color and fill by Type and DietGroup
# NOTE: DietGroup vectors are in order: c(dark, medium, light, lightest)
hex_cec <- "#2d004b" # purple - dark
hex_prx <- c("#542788", "#f2ebf9") # purple - medium, lightest
hex_dtl <- "#8073ac" # purple - light
hex_fec <- "#8c510a" # brown
hex_C <- c("#2d8c27", "#33a02c", "#9fe39b", "#ecf9eb") # green
hex_R <- c("#1a6698", "#1f78b4", "#92c8ec", "#e9f4fb") # blue
hex_F <- c("#b24b01", "#e66101", "#feb580", "#fff0e6")  # orange

# ggplot2 codes for linetypes by DietGroup
lne_C <- "dotdash"
lne_R <- "dashed"
lne_F <- "solid"

# create a vector naming all of the above (useful when saving workspaces)
SarScroll <- c("shp_cec", "shp_prx", "shp_dtl", "shp_fec",
               "shp_C", "shp_R", "shp_F",
               "hex_cec", "hex_prx", "hex_dtl", "hex_fec",
               "hex_C", "hex_R", "hex_F",
               "lne_C", "lne_R", "lne_F")

# save universal objects workspace
COSMOS <- c(cosmos, "COSMOS", SarScroll)
save(list = COSMOS, file = ofv_COSMOS_wksp)

### ************************************
### ^^^I - MACHINE/PACKAGE/VERSION Info ----
### ************************************

# ^^^ = check "# capture R package-related information:" to ensure it is correct

# NOTE: section I requires objects from the PREFACE to be in the environment
# ** note for KDP: section I version 0.2 ** #

# R packages accessed via require:
require(ggplot2, quietly = T)
require(ggpubr, quietly = T)

# R packages accessed via namespace:
# benchmarkme
# dplyr
# DECIPHER
# Biostrings
# phangorn
# ape
# ggbiplot
# GUniFrac
# vegan
# zCompositions
# PMCMR
# reshape2

# capture R package-related information:
I.Rpac_ctg <- "R package version"
I.Rpackge_a <- data.frame(info = "benchmarkme", section = "I", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_b <- data.frame(info = "dplyr", section = "I; A; S; B; E", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_c <- data.frame(info = "ape", section = "T; B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_d <- data.frame(info = "Biostrings", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_e <- data.frame(info = "DECIPHER", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_f <- data.frame(info = "phangorn", section = "T", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_g <- data.frame(info = "ggplot2", section = "B; E", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_h <- data.frame(info = "ggpubr", section = "B; E", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_i <- data.frame(info = "ggbiplot", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_j <- data.frame(info = "GUniFrac", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_k <- data.frame(info = "vegan", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_l <- data.frame(info = "zCompositions", section = "B", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_m <- data.frame(info = "PMCMR", section = "E", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_n <- data.frame(info = "reshape2", section = "E", 
                          category = I.Rpac_ctg, script = name_scrpt, 
                          stringsAsFactors = F)
I.Rpackge_0 <- rbind(I.Rpackge_a, I.Rpackge_b, I.Rpackge_c, I.Rpackge_d, 
                     I.Rpackge_e, I.Rpackge_f, I.Rpackge_g, I.Rpackge_h, 
                     I.Rpackge_i, I.Rpackge_j, I.Rpackge_k, I.Rpackge_l,
                     I.Rpackge_m, I.Rpackge_n)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: captureVersions() version 0.1 ** #
# function takes an input data.frame with:
# column 'info' listing the package name and returns a data.frame with a ...
# ... new column 'value' listing the version information for the package
# NOTE: input column names are rigid; column 'info' must list package name
captureVersions <- function(data = data.frame) {
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  # store original input df and format the new df
  new_dat <- data
  new_dat$value <- ""
  for(i in 1:length(new_dat$info)) {
    # create package name variable to be tested to ensure correct inputs/outputs
    pack <- unlist(
      packageDescription(
        new_dat$info[i], fields = c("Package", "Version"))[1])
    # test if package name in input data matches 'pack' variable
    # if TRUE, add package version value to the new data.frame
    if (pack == new_dat$info[i]) {
      new_dat$value[i] <- unlist(
        packageDescription(
          new_dat$info[i], fields = c("Package", "Version"))[2])
    }
    # if FALSE, print error message
    else {
      if (!pack == new_dat$info[i]) {
        stop("'pack' variable returned by packageDescription();
       does not match row value in column 'package' of input data;
       incorrect package version valuermation likely returned")
      }
    }
  }
  return(new_dat)
}
# 
# # example usage:
# new_dataframe <- captureVersions(data = dataframe)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# captureVersions() function
I.Rpackge_1 <- captureVersions(I.Rpackge_0)
I.Rpackge <- dplyr::select(I.Rpackge_1, category, info, value, script, section)

# capture Project-related information:
I.project <- data.frame(category = "project", 
                        info = "name of project directory",
                        value = name_projd,
                        script = name_scrpt,
                        section = "all", stringsAsFactors = F)

# capture PATH-related information:
I.path_ctg <- "filepath"
I.path_sec <- "all"
I.pathsto_a <- data.frame(info = "working directory",
                          value = paste(getwd(), "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_b <- data.frame(info = "path to project directory",
                          value = paste(path_projd, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_c <- data.frame(info = "path to project's central '/vault/'",
                          value = paste(path_vault, "/", sep = ""),
                          category = I.path_ctg, script = name_scrpt,
                          section = I.path_sec, stringsAsFactors = F)
I.pathsto_0 <- rbind(I.pathsto_a, I.pathsto_b, I.pathsto_c)
I.pathsto <- dplyr::select(I.pathsto_0, category, info, value, script, section)


# capture Machine-related information:
I.mach_ctg <- "machine"
I.mach_sec <- "all"
I.mach_cvr <- 1073741824 # number of bytes in 1GB (used for conversion of RAM)
I.machine_a <- data.frame(info = "OS",
                          value = sessionInfo()$running,
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_b <- data.frame(info = "processor",
                          value = benchmarkme::get_cpu()$model_name,
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_c <- data.frame(info = "number of cores",
                          value = parallel::detectCores(logical = F),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_d <- data.frame(info = "number of threads",
                          value = parallel::detectCores(logical = T),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_e <- data.frame(info = "RAM",
                          value = 
                            paste(
                              as.numeric(benchmarkme::get_ram()) / I.mach_cvr,
                              "GB", sep = ""),
                          category = I.mach_ctg, script = name_scrpt,
                          section = I.mach_sec, stringsAsFactors = F)
I.machine_0 <- rbind(I.machine_a, I.machine_b, I.machine_c, I.machine_d, 
                     I.machine_e)
I.machine <- dplyr::select(I.machine_0, category, info, value, script, section)

# capture R-related information:
I.rlan_ctg <- "base R"
I.rlan_sec <- "all"
I.baseRpl_a <- data.frame(info = "version",
                          value = R.Version()$version.string,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_b <- data.frame(info = "nickname",
                          value = R.Version()$nickname,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_c <- data.frame(info = "platform",
                          value = R.Version()$platform,
                          category = I.rlan_ctg, script = name_scrpt,
                          section = I.rlan_sec, stringsAsFactors = F)
I.baseRpl_0 <- rbind(I.baseRpl_a, I.baseRpl_b, I.baseRpl_c)
I.baseRpl <- dplyr::select(I.baseRpl_0, category, info, value, script, section)

# capture RStudio-related information:
I.RStudio <- data.frame(info = "version",
                        value = as.character(RStudio.Version()$version),
                        category = "R Studio", script = name_scrpt,
                        section = "all", stringsAsFactors = F)

# rbind all of the above data.frames together and write outputs
I.info <- rbind(I.project, I.machine, I.pathsto, I.baseRpl, I.RStudio, 
                I.Rpackge)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = I.info, file = I.ofv_info)

I.obj <- ls(pattern = "I.")
I.lst <- c(I.obj[grep(pattern = "I.", x = I.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, "captureVersions")
save(list = I.lst, file = I.ofv_wksp)

### ************************************
### ^^^DESCRIPTIONS FOR EACH Section ----
### ************************************

# ^^^ descriptions not complete for sections outside of this script

## SECTION A - create a master table (EBTKS)
# purpose: take outputs from QIIME 2 and combine them into a master table...
# ... to serve as the entry point for all downstream processing/analysis
# an EBTKS feature table contains ...
# ... hashed FeatureIDs & representative sequences for all ASVs ...
# ... full & truncated Greengenes & SILVA taxonomic lineages for all ASVs ...
# ... absolute counts for all ASVs in each sample of the dataset
# NOTE: EBTKS = Everything But The Kitchen Sink

## SECTION S - subset the processed EBTKS table to isolate samples of interest
# purpose: subset tables with specific groups are needed for a number of ...
# .. downstream analyses

## SECTION T - create a rooted phylogenetic tree from processed EBTKS table
# purpose: a rooted tree is needed for phylogenetic-based analyses
# NOTE: be aware that this is a [typically] computationally intensive section

## Section B - beta diversity and microbiota composition
#
## Ssection C - colonization (i.e. ASV/taxon relative abundances)
#
## Section D - differential abundance testing
# purpose:
# this section tests for differentially abundant features (ASVs) b/w two groups
## Section L - tumor outcomes - neoplastic lesions
#

### ************************************
### A - FUNCTIONS ----
### ************************************

# create a vector naming all Section A functions (used when saving workspace)
A.functions <- c("fasta_to_df", "trunc_tax", "testif_fts_lost", 
                 "gather_lost_fts")

# ** note for KDP: fasta_to_df() version 0.2 ** #
# function takes an input .fasta file in two-line format which looks like:
# line 1: >header
# line 2: DNA sequence
# ... and creates a two column dataframe
# fasta_to_df == fasta file to data.frame
fasta_to_df <- function(fasta_file, hdr_colname = "", seq_colname = "") {
  # read in the .fasta file line by line
  fasta <- readLines(fasta_file)
  # identify header lines by finding line numbers with the '>' character
  hlines <- grep(pattern = ">", x = fasta)
  # create a three column df containing:
  # header line numbers (hdr)
  # sequence line number beginning (beg)
  # sequence line number end (end)
  slines <- data.frame(hdr = hlines,
                       beg = hlines+1,
                       end = c((hlines-1)[-1], length(fasta)))
  # create a vector of identical length to hdr_lines to be used for
  # storing values obtained in the loop below
  storage_vec <- rep(NA, length(hlines))
  # loop to obtain sequences
  for(i in 1:length(hlines)) {
    storage_vec[i] <- paste(fasta[slines$beg[i]:slines$end[i]], collapse = "")
  }
  # create a two column df containing:
  # the header with '>' character replaced (hdr)
  # representative sequence associated with the header (seq)
  new_df_v0 <- data.frame("hdr" = gsub(pattern = ">", replacement = "",
                                       x = fasta[hlines]),
                          "seq" = storage_vec)
  # replace generic column names with names specified in the function input
  new_df_v1 <- new_df_v0
  hdr_colnum <- which(names(new_df_v1) == "hdr")
  seq_colnum <- which(names(new_df_v1) == "seq")
  names(new_df_v1)[hdr_colnum] <- hdr_colname
  names(new_df_v1)[seq_colnum] <- seq_colname
  # cols are currently factors, coerce them to characters
  new_df_v2 <- new_df_v1
  new_df_v2[, hdr_colnum] <- as.character(new_df_v2[, hdr_colnum])
  new_df_v2[, seq_colnum] <- as.character(new_df_v2[, seq_colnum])
  return(new_df_v2)
}
# 
# # example usage:
# df_fasta <- fasta_to_df(fasta_file = ifp_fasta,
#                         hdr_colname = "FeatureID",
#                         seq_colname = "RepSeq")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: trunc_tax() version 0.4 ** #
# trunc_tax() function requires an input df with a column containing either:
# greengenes formatted taxonomic strings, e.g.:
# k__; p__ ; c__; o__; f__; g__; s__"
# or silva taxonomic strings, e.g.:
# "D_0__;D_1__;D_2__;D_3__;D_4__;D_5__;D_6__"
# the output df will contain a new column with:
# the desired taxonomic level for the specified database naming convention
# NOTE: if the above strings are not in consecutive order...
# ... the output may not be correct ...
# i.e. if string input is missing D_2__ = "D_0__;D_1__Example;D_3__Example;D_4__"
# output for ;D_1__ = "Example;D_3__Example;D_4__" rather than just "Example"
# output for ;D_2__ = "Unassigned"
# output for ;D_3__ = "Example"
# trunc_tax = truncate taxonomic string
trunc_tax <- function(data = data.frame, tax_type = c("Greengenes", "SILVA"),
                      tax_lvl = c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species"),
                      input_col = "", str1 = "", str2 = "", classifier = "") {
  
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA") {
    stop("tax_type should be either 'Greengenes' or 'SILVA'")
  }
  if (!tax_lvl == "Kingdom" && !tax_lvl == "Phylum" && !tax_lvl == "Class" &&
      !tax_lvl == "Order" && !tax_lvl == "Family" && !tax_lvl == "Genus" &&
      !tax_lvl == "Species") {
    stop("tax_lvl should be one of: 'Kingdom', 'Phylum', 'Class', 'Order',
        'Family', 'Genus', 'Species'")
  }
  if (!inherits(input_col, "character")) {
    stop("input_col must be character")
  }
  if (!inherits(str1, "character")) {
    stop("str1 must be character")
  }
  if (!inherits(str2, "character")) {
    stop("str2 must be character")
  }
  if (!inherits(classifier, "character")) {
    stop("classifier must be character")
  }
  
  # if tax_lvl is Species, pasting the string split is handled differently
  num <- ifelse(tax_lvl == "Species", yes = 2, no = 1)
  
  # store original input df and format the new df
  new_df1 <- data
  new_df1$new_col <- 1
  
  # loop through data and truncate the taxonomic strings
  for(row in 1:nrow(new_df1)) {
    # split the string in the specified input column with input of str1
    splt1 <- base::strsplit(new_df1[row, input_col], str1)[[1]][2]
    # goal of first ifelse test:
    # yes = tax level was not assigned
    # no = tax level was likely assigned, split the string with input of str2
    splt2 <- ifelse(is.na(splt1), yes = "Unassigned",
                    no = paste("", base::strsplit(splt1, str2)[[1]][num], 
                               sep = ""))
    # goal of second ifelse test:
    # yes = tax level was not assigned
    # no = tax level was actually assigned
    new_df1[row, "new_col"] <- ifelse(nchar(splt2) == 0 | "NA" %in% splt2,
                                      yes = "Unassigned", no = splt2)
    # if tax_type is Greengenes, tax_lvl is Species and row is not Unassigned,
    # append Genus in front of Species assignment
    if (tax_type == "Greengenes" & tax_lvl == "Species" &
        !new_df1[row, "new_col"] == "Unassigned") {
      new_df1[row, "new_col"] <- paste(base::strsplit(splt1, str2)[[1]][1],
                                       base::strsplit(splt1, str2)[[1]][2], 
                                       sep = " ")
    }
  } # close loop
  
  # internal warning in the event that all rows in the new column are Unassigned
  # this outcome could be correct if the specified tax lvl had no assignemnts;
  # however, if assignments are expected for the specified tax lvl, this outcome
  # is incorrect and the function inputs need to be double checked
  if (isTRUE(all(grepl(pattern = "Unassigned", x = new_df1[, "new_col"],
                       ignore.case = FALSE, fixed = TRUE)))) {
    warning("in new_col all row values in output column are 'Unassigned'")
  }
  
  # format the new df
  new_df2 <- new_df1
  new_col_num <- which(names(new_df2) == 'new_col')
  if (tax_type == "Greengenes") {
    names(new_df2)[new_col_num] <- paste(classifier, "ggs", tax_lvl, sep = ".")
  }
  if (tax_type == "SILVA") {
    names(new_df2)[new_col_num] <- paste(classifier, "slv", tax_lvl, sep = ".")
  }
  return(new_df2)
}
#
# example usage:
# inputs are incorrect, prints warning
#new_df <- trunc_tax(data = a.df, tax_type = "silva",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

# corrected inputs, no warning
#new_df <- trunc_tax(data = a.df, tax_type = "greengenes",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.1 ** #
# testif_fts_lost() and gather_lost_fts() functions take identical inputs:
# 1: a raw taxonomy data.frame with column names changed (see: A - STEP 1a)
# 2: a processed taxanomy data frame from one of A - Step 1a; 1b; 1c
# 3: the type specificaion for the input to test_df
# 4: the taxonomic database used for the input to test_df
# 5: input column name to join by (e.g. FeatureID or some equivalent)
# NOTE: input data.frames must be from the same taxonomy database
# NOTE: I understand that there are a million ways to do this, I chose one.
# NOTE: both functions are dependent upon package dplyr()
# ** note for KDP: future versions to halt incorrect _type combos ** #

# the testif_fts_lost() function: 
# checks if any features were lost after a processing step was performed
# if features were lost, an error message is printed and ...
# ... an empty data.frame is returned in order to prevent crucial parts ...
# of the downstream code from continuing on
# the gather_lost_fts() function: 
# identifies any lost features after a processing step was performed
# typically, the SILVA database is the issue with features missing when the ...
# _slv_lws variable is created from the _slv_trunc_fix variable in:
# A - STEP 1c - format taxonomy tables: determine lowest assignment level
# this specific error can be overcome by adding to the _slv_trunc_fix variable
# # more details about that process can be found in:
# A - STEP 1b - format taxonomy tables: "fix" truncated lineages

testif_fts_lost <- function(raw_df = data.frame, test_df = data.frame, 
                            test_type = c("trunc", "trunc_fix", "lws", 
                                          "frm", "mrg"),
                            tax_type = c("Greengenes", "SILVA", "merging"),
                            col_join = "") {
  # internal checks to ensure correct input classes
  if (!inherits(raw_df, "data.frame")) {
    stop("input raw_df must be class 'data.frame'")
  }
  if (!inherits(test_df, "data.frame")) {
    stop("input test_df must be class 'data.frame'")
  }
  if (!test_type == "trunc" && !test_type == "trunc_fix" && 
      !test_type == "lws" && !test_type == "frm" && !test_type == "mrg") {
    stop("test_type must be one of: 'trunc', 'trunc_fix', 'lws', 'frm', 'mrg'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA" && 
      !tax_type == "merging") {
    stop("tax_type must be one of: 'Greengenes' or 'SILVA', 'merging'")
  }
  if (!inherits(col_join, "character")) {
    stop("col_join must be character")
  }
  
  # create the error message for !test_type == "mrg" and !tax_type == "merging"
  err_pr1 <- "\ndetected lost fts in '_"
  err_pr2 <- "' data.frame from "
  err <- paste(err_pr1, test_type, err_pr2, tax_type, sep = "")
  
  # create the error message for test_type == "mrg" and tax_type == "merging"
  err_prg <- "' data.frame after database "
  err_mrg <- paste(err_pr1, test_type, err_prg, tax_type, sep = "")
  
  # add the final part to the error message
  err_fpr <- "run: gather_lost_fts() function to identify missing features"
  fin_err <- paste(err, err_fpr, sep = "\n")
  fin_err_mrg <- paste(err_mrg, err_fpr, sep = "\n")
  
  # anti_join the two input data.frames to determine if features were lost
  anti_df <- dplyr::anti_join(x = raw_df, y = test_df, by = col_join)
  
  # using an ifelse statement, create a logical variable to test the anti_df
  LOGICAL_anti_df <- ifelse(nrow(anti_df) == 0, yes = T, no = F)
  
  # test the LOGICAL variable; 
  # if FALSE, return an empty data.frame which will halt downstream code and ...
  # ... print the appropriate error message based upon the input to test_type
  
  if (!test_type == "merging") {  
    if (!isTRUE(LOGICAL_anti_df)) {
      test_df <- data.frame()
      stop(fin_err)
    }
  }
  
  if (test_type == "merging") {  
    if (!isTRUE(LOGICAL_anti_df)) {
      test_df <- data.frame()
      stop(fin_err_mrg)
    }
  }
  return(test_df)
}

gather_lost_fts <- function(raw_df = data.frame, test_df = data.frame, 
                            test_type = c("trunc", "trunc_fix", "lws", 
                                          "frm", "mrg"),
                            tax_type = c("Greengenes", "SILVA", "merging"),
                            col_join = "") {
  # internal checks to ensure correct input classes
  if (!inherits(raw_df, "data.frame")) {
    stop("input raw_df must be class 'data.frame'")
  }
  if (!inherits(test_df, "data.frame")) {
    stop("input test_df must be class 'data.frame'")
  }
  if (!test_type == "trunc" && !test_type == "trunc_fix" && 
      !test_type == "lws" && !test_type == "frm" && !test_type == "mrg") {
    stop("test_type must be one of: 'trunc', 'trunc_fix', 'lws', 'frm', 'mrg'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA" && 
      !tax_type == "merging") {
    stop("tax_type must be one of: 'Greengenes' or 'SILVA', 'merging'")
  }
  if (!inherits(col_join, "character")) {
    stop("col_join must be character")
  }
  
  # anti_join the two input data.frames to determine if features were lost
  anti_df <- dplyr::anti_join(x = raw_df, y = test_df, by = col_join)
  return(anti_df)
}
#
# example usage:
# A.tst_slv_lws <- testif_fts_lost(raw_df = A.int_slv, test_df = A.int_slv_lws, 
#                                   test_type = "lws", tax_type = "SILVA",
#                                   col_join = "FeatureID")
# A.lost_slv_lws <- gather_lost_fts(raw_df = A.int_slv, test_df = A.int_slv_lws, 
#                                   test_type = "lws", tax_type = "SILVA",
#                                   col_join = "FeatureID")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### A - STEP 1a - format taxonomy tables: truncate lineages ----
### ************************************

# ** note for KDP: AS truncate lineages version 0.4 ** #

# this step truncates Greengenes & SILVA taxonomic lineages, placing them...
# ...into new cols: .Kingdom; .Phylum; .Class; .Order; .Family; .Genus; .Species
# this step also determines the lowest assigned level, placing that info...
# ... into new cols: .lws.txn and .lws.lvl

# read in taxonomy tables
A.raw_int_ggs <- read.table(A.ifv_int_ggs, header = T, sep = "\t", as.is = T,
                            stringsAsFactors = F)
A.raw_int_slv <- read.table(A.ifv_int_slv, header = T, sep = "\t", as.is = T,
                            stringsAsFactors = F, comment.char = "", quote = "")

# rename columns
A.int_ggs <- dplyr::rename(A.raw_int_ggs, FeatureID = Feature.ID,
                           int.ggs.tax = Taxon, int.ggs.cnf = Confidence)
A.int_slv <- dplyr::rename(A.raw_int_slv, FeatureID = Feature.ID,
                           int.slv.tax = Taxon, int.slv.cnf = Confidence)

# truncate taxonomic lineages for each level, and reorder columns
A.int_ggs_L1 <- trunc_tax(data = A.int_ggs, tax_type = "Greengenes",
                          tax_lvl = "Kingdom", input_col = "int.ggs.tax",
                          str1 = "k__", str2 = "; p__", classifier = "int")
A.int_ggs_L2 <- trunc_tax(data = A.int_ggs_L1, tax_type = "Greengenes",
                          tax_lvl = "Phylum", input_col = "int.ggs.tax",
                          str1 = "; p__", str2 = "; c__", classifier = "int")
A.int_ggs_L3 <- trunc_tax(data = A.int_ggs_L2, tax_type = "Greengenes",
                          tax_lvl = "Class", input_col = "int.ggs.tax",
                          str1 = "; c__", str2 = "; o__", classifier = "int")
A.int_ggs_L4 <- trunc_tax(data = A.int_ggs_L3, tax_type = "Greengenes",
                          tax_lvl = "Order", input_col = "int.ggs.tax",
                          str1 = "; o__", str2 = "; f__", classifier = "int")
A.int_ggs_L5 <- trunc_tax(data = A.int_ggs_L4, tax_type = "Greengenes",
                          tax_lvl = "Family", input_col = "int.ggs.tax",
                          str1 = "; f__", str2 = "; g__", classifier = "int")
A.int_ggs_L6 <- trunc_tax(data = A.int_ggs_L5, tax_type = "Greengenes",
                          tax_lvl = "Genus", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A.int_ggs_L7 <- trunc_tax(data = A.int_ggs_L6, tax_type = "Greengenes",
                          tax_lvl = "Species", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A.int_ggs_trunc <- dplyr::select(A.int_ggs_L7, FeatureID, int.ggs.tax,
                                 int.ggs.Kingdom, int.ggs.Phylum, int.ggs.Class,
                                 int.ggs.Order, int.ggs.Family, int.ggs.Genus,
                                 int.ggs.Species, int.ggs.cnf)

A.int_slv_L1 <- trunc_tax(data = A.int_slv, tax_type = "SILVA",
                          tax_lvl = "Kingdom", input_col = "int.slv.tax",
                          str1 = "D_0__", str2 = ";D_1__", classifier = "int")
A.int_slv_L2 <- trunc_tax(data = A.int_slv_L1, tax_type = "SILVA",
                          tax_lvl = "Phylum", input_col = "int.slv.tax",
                          str1 = ";D_1__", str2 = ";D_2__", classifier = "int")
A.int_slv_L3 <- trunc_tax(data = A.int_slv_L2, tax_type = "SILVA",
                          tax_lvl = "Class", input_col = "int.slv.tax",
                          str1 = ";D_2__", str2 = ";D_3__", classifier = "int")
A.int_slv_L4 <- trunc_tax(data = A.int_slv_L3, tax_type = "SILVA",
                          tax_lvl = "Order", input_col = "int.slv.tax",
                          str1 = ";D_3__", str2 = ";D_4__", classifier = "int")
A.int_slv_L5 <- trunc_tax(data = A.int_slv_L4, tax_type = "SILVA",
                          tax_lvl = "Family", input_col = "int.slv.tax",
                          str1 = ";D_4__", str2 = ";D_5__", classifier = "int")
A.int_slv_L6 <- trunc_tax(data = A.int_slv_L5, tax_type = "SILVA",
                          tax_lvl = "Genus", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A.int_slv_L7 <- trunc_tax(data = A.int_slv_L6, tax_type = "SILVA",
                          tax_lvl = "Species", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A.int_slv_trunc <- dplyr::select(A.int_slv_L7, FeatureID, int.slv.tax,
                                 int.slv.Kingdom, int.slv.Phylum, int.slv.Class,
                                 int.slv.Order, int.slv.Family, int.slv.Genus,
                                 int.slv.Species, int.slv.cnf)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in either of the truncated lineage dfs

A.test_ggs_trunc <- testif_fts_lost(raw_df = A.int_ggs, 
                                    test_df = A.int_ggs_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "Greengenes",
                                    col_join = "FeatureID")
A.test_slv_trunc <- testif_fts_lost(raw_df = A.int_slv, 
                                    test_df = A.int_slv_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "SILVA",
                                    col_join = "FeatureID")
A.lost_ggs_trunc <- gather_lost_fts(raw_df = A.int_ggs, 
                                    test_df = A.int_ggs_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "Greengenes",
                                    col_join = "FeatureID")
A.lost_slv_trunc <- gather_lost_fts(raw_df = A.int_slv, 
                                    test_df = A.int_slv_trunc, 
                                    test_type = "trunc", 
                                    tax_type = "SILVA",
                                    col_join = "FeatureID")

### ************************************
### A - STEP 1b - format taxonomy tables: "fix" truncated lineages ----
### ************************************

# ** note for KDP: AS "fix" truncated lineages version 0.9 ** #

# An aside on the subject of truncated taxonomic lineages:

# truncated taxonomic lineages are great and mostly useful; however...
# ... a variety of SILVA assignments lose useful information when truncated
# e.g. Order = Lineage IV
# e.g. Family = uncultured Rubrobacteria bacterium
# e.g. Genus = AAP99
# this is especially problematic for the SILVA assignment Family = Family XI ...
# ... this family may be in the Order Bacillales or the Order Clostridiales
# the goal here is to add info to make truncated assignments useful/meaningful
# e.g. Family XI has info prepended for the Order
# e.g. Order: Lineage IV is changed to Elusimicrobia Lineage IV
# NOTE: the primary focus here is on SILVA as this database is actually updated
# NOTE: these are "running-lists" and are likely dataset specific...
# ... i.e. not all lineages will always be present in every dataset
# NOTE: dtabase version info is located in COMPUTER/VERSION/PACKAGE INFO

# end of aside, begin comments for the code below

# create vectors to be used to index column names
# NOTE: these vectors will be indexed using numbers; DO NOT ALTER THEIR ORDER
A.ggs_col <- c("int.ggs.Kingdom", "int.ggs.Phylum", "int.ggs.Class",
               "int.ggs.Order", "int.ggs.Family","int.ggs.Genus",
               "int.ggs.Species")

A.slv_col <- c("int.slv.Kingdom", "int.slv.Phylum", "int.slv.Class",
               "int.slv.Order", "int.slv.Family","int.slv.Genus",
               "int.slv.Species")

# the below is performed on one assignment at a time; the process is as follows:
# (1) create a vector to filter the lineage that needs "fixed"
# (2) isolate that lineage from the A.test_slv_trunc df created in A - STEP 1A
# (3) fix the lineage
# (4) rinse and repeat for all lineages that need to be "fixed"
# NOTE: the capital letter combinations indicate taxonomic level...
# ... i.e. PC = Phylum-Class; OFG = Order-Family-Genus, etc.

# single level fixes first:
# Phylum-Class
A.slv_vec_PC_01 <- c("Acidobacteria", "Subgroup 6")
A.slv_flt_PC_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[2]]] == A.slv_vec_PC_01[1],
                                 .data[[A.slv_col[3]]] == A.slv_vec_PC_01[2])
A.slv_fix_PC_01 <- dplyr::mutate(
  A.slv_flt_PC_01,
  int.slv.Class = replace(x = .data[[A.slv_col[2]]],
                          list = .data[[A.slv_col[3]]] == A.slv_vec_PC_01[1] &
                            .data[[A.slv_col[2]]] == A.slv_vec_PC_01[2],
                          values = paste(A.slv_vec_PC_01[1],
                                         A.slv_vec_PC_01[2],
                                         sep = " ")))

# Class-Order
A.slv_vec_CO_01 <- c("Anaerolineae", "SBR1031")
A.slv_flt_CO_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_01[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_01[2])
A.slv_fix_CO_01 <- dplyr::mutate(
  A.slv_flt_CO_01,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_01[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_01[2],
                          values = paste(A.slv_vec_CO_01[1],
                                         A.slv_vec_CO_01[2],
                                         sep = " ")))

A.slv_vec_CO_02 <- c("Dehalococcoidia", "S085")
A.slv_flt_CO_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_02[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_02[2])
A.slv_fix_CO_02 <- dplyr::mutate(
  A.slv_flt_CO_02,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_02[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_02[2],
                          values = paste(A.slv_vec_CO_02[1],
                                         A.slv_vec_CO_02[2],
                                         sep = " ")))

A.slv_vec_CO_03 <- c("Elusimicrobia", "Lineage IV")
A.slv_flt_CO_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[3]]] == A.slv_vec_CO_03[1],
                                 .data[[A.slv_col[4]]] == A.slv_vec_CO_03[2])
A.slv_fix_CO_03 <- dplyr::mutate(
  A.slv_flt_CO_03,
  int.slv.Class = replace(x = .data[[A.slv_col[3]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_CO_03[1] &
                            .data[[A.slv_col[4]]] == A.slv_vec_CO_03[2],
                          values = paste(A.slv_vec_CO_03[1],
                                         A.slv_vec_CO_03[2],
                                         sep = " ")))

# Order-Family
A.slv_vec_OF_01 <- c("Bacillales", "Family XI")
A.slv_flt_OF_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_01[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_01[2])
A.slv_fix_OF_01 <- dplyr::mutate(
  A.slv_flt_OF_01,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_01[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_01[2],
                           values = paste(A.slv_vec_OF_01[1],
                                          A.slv_vec_OF_01[2],
                                          sep = " ")))

A.slv_vec_OF_02 <- c("Bacteroidales", "p-2534-18B5 gut group")
A.slv_flt_OF_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_02[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_02[2])
A.slv_fix_OF_02 <- dplyr::mutate(
  A.slv_flt_OF_02,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_02[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_02[2],
                           values = paste(A.slv_vec_OF_02[1],
                                          A.slv_vec_OF_02[2],
                                          sep = " ")))

A.slv_vec_OF_03 <- c("Clostridiales", "Family XI")
A.slv_flt_OF_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_03[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_03[2])
A.slv_fix_OF_03 <- dplyr::mutate(
  A.slv_flt_OF_03,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_03[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_03[2],
                           values = paste(A.slv_vec_OF_03[1],
                                          A.slv_vec_OF_03[2],
                                          sep = " ")))

A.slv_vec_OF_04 <- c("Flavobacteriales", "NS9 marine group")
A.slv_flt_OF_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_04[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_04[2])
A.slv_fix_OF_04 <- dplyr::mutate(
  A.slv_flt_OF_04,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_04[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_04[2],
                           values = paste(A.slv_vec_OF_04[1],
                                          A.slv_vec_OF_04[2],
                                          sep = " ")))

A.slv_vec_OF_05 <- c("Kryptoniales", "BSV26")
A.slv_flt_OF_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_05[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_05[2])
A.slv_fix_OF_05 <- dplyr::mutate(
  A.slv_flt_OF_05,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_05[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_05[2],
                           values = paste(A.slv_vec_OF_05[1],
                                          A.slv_vec_OF_05[2],
                                          sep = " ")))

A.slv_vec_OF_06 <- c("Myxococcales", "BIrii41")
A.slv_flt_OF_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_06[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_06[2])
A.slv_fix_OF_06 <- dplyr::mutate(
  A.slv_flt_OF_06,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_06[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_06[2],
                           values = paste(A.slv_vec_OF_06[1],
                                          A.slv_vec_OF_06[2],
                                          sep = " ")))

A.slv_vec_OF_07 <- c("Solirubrobacterales", "67-14")
A.slv_flt_OF_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_07[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_07[2])
A.slv_fix_OF_07 <- dplyr::mutate(
  A.slv_flt_OF_07,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_07[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_07[2],
                           values = paste(A.slv_vec_OF_07[1],
                                          A.slv_vec_OF_07[2],
                                          sep = " ")))

A.slv_vec_OF_08 <- c("Thermoanaerobacterales", "Family III")
A.slv_flt_OF_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_08[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_08[2])
A.slv_fix_OF_08 <- dplyr::mutate(
  A.slv_flt_OF_08,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_08[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_08[2],
                           values = paste(A.slv_vec_OF_08[1],
                                          A.slv_vec_OF_08[2],
                                          sep = " ")))

A.slv_vec_OF_09 <- c("Thermomicrobiales", "JG30-KF-CM45")
A.slv_flt_OF_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_09[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_09[2])
A.slv_fix_OF_09 <- dplyr::mutate(
  A.slv_flt_OF_09,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_09[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_09[2],
                           values = paste(A.slv_vec_OF_09[1],
                                          A.slv_vec_OF_09[2],
                                          sep = " ")))

A.slv_vec_OF_10 <- c("Clostridiales", "Family XIII")
A.slv_flt_OF_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_10[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_10[2],
                                 !grepl(pattern = A.slv_vec_OF_10[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_10 <- dplyr::mutate(
  A.slv_flt_OF_10,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_10[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_10[2],
                           values = paste(A.slv_vec_OF_10[1],
                                          A.slv_vec_OF_10[2],
                                          sep = " ")))

A.slv_vec_OF_11 <- c("Verrucomicrobiales", "DEV007")
A.slv_flt_OF_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_11[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_11[2],
                                 !grepl(pattern = A.slv_vec_OF_11[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_11 <- dplyr::mutate(
  A.slv_flt_OF_11,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_11[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_11[2],
                           values = paste(A.slv_vec_OF_11[1],
                                          A.slv_vec_OF_11[2],
                                          sep = " ")))

A.slv_vec_OF_12 <- c("Myxococcales", "P3OB-42")
A.slv_flt_OF_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[4]]] == A.slv_vec_OF_12[1],
                                 .data[[A.slv_col[5]]] == A.slv_vec_OF_12[2],
                                 !grepl(pattern = A.slv_vec_OF_12[2],
                                        x = .data[[A.slv_col[6]]]))
A.slv_fix_OF_12 <- dplyr::mutate(
  A.slv_flt_OF_12,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OF_12[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OF_12[2],
                           values = paste(A.slv_vec_OF_12[1],
                                          A.slv_vec_OF_12[2],
                                          sep = " ")))

# Family-Genus
A.slv_vec_FG_01 <- c("Actinomycetaceae", "F0332")
A.slv_flt_FG_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_01[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_01[2])
A.slv_fix_FG_01 <- dplyr::mutate(
  A.slv_flt_FG_01,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_01[2],
                          values = paste(A.slv_vec_FG_01[1],
                                         A.slv_vec_FG_01[2],
                                         sep = " ")))

A.slv_vec_FG_02 <- c("Burkholderiaceae", "AAP99")
A.slv_flt_FG_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_02[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_02[2])
A.slv_fix_FG_02 <- dplyr::mutate(
  A.slv_flt_FG_02,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_02[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_02[2],
                          values = paste(A.slv_vec_FG_02[1],
                                         A.slv_vec_FG_02[2],
                                         sep = " ")))

A.slv_vec_FG_03 <- c("Lachnospiraceae", "A2")
A.slv_flt_FG_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_03[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_03[2])
A.slv_fix_FG_03 <- dplyr::mutate(
  A.slv_flt_FG_03,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_03[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_03[2],
                          values = paste(A.slv_vec_FG_03[1],
                                         A.slv_vec_FG_03[2],
                                         sep = " ")))

A.slv_vec_FG_04 <- c("Lachnospiraceae", "ASF356")
A.slv_flt_FG_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_04[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_04[2])
A.slv_fix_FG_04 <- dplyr::mutate(
  A.slv_flt_FG_04,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_04[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_04[2],
                          values = paste(A.slv_vec_FG_04[1],
                                         A.slv_vec_FG_04[2],
                                         sep = " ")))

A.slv_vec_FG_05 <- c("Lachnospiraceae", "CAG-56")
A.slv_flt_FG_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_05[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_05[2])
A.slv_fix_FG_05 <- dplyr::mutate(
  A.slv_flt_FG_05,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_05[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_05[2],
                          values = paste(A.slv_vec_FG_05[1],
                                         A.slv_vec_FG_05[2],
                                         sep = " ")))

A.slv_vec_FG_06 <- c("Lachnospiraceae", "GCA-900066575")
A.slv_flt_FG_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_06[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_06[2])
A.slv_fix_FG_06 <- dplyr::mutate(
  A.slv_flt_FG_06,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_06[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_06[2],
                          values = paste(A.slv_vec_FG_06[1],
                                         A.slv_vec_FG_06[2],
                                         sep = " ")))

A.slv_vec_FG_07 <- c("Lachnospiraceae", "GCA-900066755")
A.slv_flt_FG_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_07[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_07[2])
A.slv_fix_FG_07 <- dplyr::mutate(
  A.slv_flt_FG_07,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_07[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_07[2],
                          values = paste(A.slv_vec_FG_07[1],
                                         A.slv_vec_FG_07[2],
                                         sep = " ")))

A.slv_vec_FG_08 <- c("Lachnospiraceae", "UC5-1-2E3")
A.slv_flt_FG_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_08[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_08[2])
A.slv_fix_FG_08 <- dplyr::mutate(
  A.slv_flt_FG_08,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_08[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_08[2],
                          values = paste(A.slv_vec_FG_08[1],
                                         A.slv_vec_FG_08[2],
                                         sep = " ")))

A.slv_vec_FG_09 <- c("Pseudohongiellaceae", "BIyi10")
A.slv_flt_FG_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_09[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_09[2])
A.slv_fix_FG_09 <- dplyr::mutate(
  A.slv_flt_FG_09,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_09[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_09[2],
                          values = paste(A.slv_vec_FG_09[1],
                                         A.slv_vec_FG_09[2],
                                         sep = " ")))

A.slv_vec_FG_10 <- c("Ruminococcaceae", "DTU089")
A.slv_flt_FG_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_10[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_10[2])
A.slv_fix_FG_10 <- dplyr::mutate(
  A.slv_flt_FG_10,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_10[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_10[2],
                          values = paste(A.slv_vec_FG_10[1],
                                         A.slv_vec_FG_10[2],
                                         sep = " ")))

A.slv_vec_FG_11 <- c("Ruminococcaceae", "GCA-900066225")
A.slv_flt_FG_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_11[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_11[2])
A.slv_fix_FG_11 <- dplyr::mutate(
  A.slv_flt_FG_11,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_11[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_11[2],
                          values = paste(A.slv_vec_FG_11[1],
                                         A.slv_vec_FG_11[2],
                                         sep = " ")))

A.slv_vec_FG_12 <- c("Ruminococcaceae", "UBA1819")
A.slv_flt_FG_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_12[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_12[2])
A.slv_fix_FG_12 <- dplyr::mutate(
  A.slv_flt_FG_12,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_12[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_12[2],
                          values = paste(A.slv_vec_FG_12[1],
                                         A.slv_vec_FG_12[2],
                                         sep = " ")))

A.slv_vec_FG_13 <- c("Thermoanaerobaculaceae", "Subgroup 10")
A.slv_flt_FG_13 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_13[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_13[2])
A.slv_fix_FG_13 <- dplyr::mutate(
  A.slv_flt_FG_13,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_13[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_13[2],
                          values = paste(A.slv_vec_FG_13[1],
                                         A.slv_vec_FG_13[2],
                                         sep = " ")))

A.slv_vec_FG_14 <- c("Thermoanaerobaculaceae", "Subgroup 23")
A.slv_flt_FG_14 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_14[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_14[2])
A.slv_fix_FG_14 <- dplyr::mutate(
  A.slv_flt_FG_14,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_14[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_14[2],
                          values = paste(A.slv_vec_FG_14[1],
                                         A.slv_vec_FG_14[2],
                                         sep = " ")))

A.slv_vec_FG_15 <- c("Xanthomonadaceae", "SN8")
A.slv_flt_FG_15 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_15[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_15[2])
A.slv_fix_FG_15 <- dplyr::mutate(
  A.slv_flt_FG_15,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_15[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_15[2],
                          values = paste(A.slv_vec_FG_15[1],
                                         A.slv_vec_FG_15[2],
                                         sep = " ")))

A.slv_vec_FG_16 <- c("Ruminococcaceae", "CAG-352")
A.slv_flt_FG_16 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_FG_16[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_FG_16[2])
A.slv_fix_FG_16 <- dplyr::mutate(
  A.slv_flt_FG_16,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FG_16[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FG_16[2],
                          values = paste(A.slv_vec_FG_16[1],
                                         A.slv_vec_FG_16[2],
                                         sep = " ")))

# Genus-Species
A.slv_vec_GS_01 <- c("Bacteroides", "bacterium NLAE-zl-H46")
A.slv_flt_GS_01 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_01[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_01[2])
A.slv_fix_GS_01 <- dplyr::mutate(
  A.slv_flt_GS_01,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_01[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_01[2],
                            values = paste(A.slv_vec_GS_01[1],
                                           A.slv_vec_GS_01[2],
                                           sep = " ")))

A.slv_vec_GS_02 <- c("Christensenellaceae R-7 group", "bacterium YE57")
A.slv_flt_GS_02 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_02[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_02[2])
A.slv_fix_GS_02 <- dplyr::mutate(
  A.slv_flt_GS_02,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_02[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_02[2],
                            values = paste(A.slv_vec_GS_02[1],
                                           A.slv_vec_GS_02[2],
                                           sep = " ")))

A.slv_vec_GS_03 <- c("Desulfovibrio", "bacterium New Zealand D")
A.slv_flt_GS_03 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_03[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_03[2])
A.slv_fix_GS_03 <- dplyr::mutate(
  A.slv_flt_GS_03,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_03[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_03[2],
                            values = paste(A.slv_vec_GS_03[1],
                                           A.slv_vec_GS_03[2],
                                           sep = " ")))

A.slv_vec_GS_04 <- c("Gaiella", "actinobacterium WWH12")
A.slv_flt_GS_04 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_04[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_04[2])
A.slv_fix_GS_04 <- dplyr::mutate(
  A.slv_flt_GS_04,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_04[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_04[2],
                            values = paste(A.slv_vec_GS_04[1],
                                           A.slv_vec_GS_04[2],
                                           sep = " ")))

A.slv_vec_GS_05 <- c("Ruminiclostridium 9",
                     "bacterium enrichment culture clone M244")
A.slv_flt_GS_05 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_05[1],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_05[2])
A.slv_fix_GS_05 <- dplyr::mutate(
  A.slv_flt_GS_05,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[6]]] == A.slv_vec_GS_05[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_05[2],
                            values = paste(A.slv_vec_GS_05[1],
                                           A.slv_vec_GS_05[2],
                                           sep = " ")))

A.slv_vec_GS_06 <- c("Lachnospiraceae", "uncultured",
                     "Frisingicoccus caecimuris", "Frisingicoccus")
A.slv_flt_GS_06 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_06[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_06[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_06[3])
A.slv_fix_GS_06 <- dplyr::mutate(
  A.slv_flt_GS_06,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_06[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_06[2],
                          values = A.slv_vec_GS_06[4]))

A.slv_vec_GS_07 <- c("Ruminococcaceae", "uncultured",
                     "Clostridium phoceensis", "Clostridium")
A.slv_flt_GS_07 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_07[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_07[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_07[3])
A.slv_fix_GS_07 <- dplyr::mutate(
  A.slv_flt_GS_07,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_07[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_07[2],
                          values = A.slv_vec_GS_07[4]))

A.slv_vec_GS_08 <- c("Lachnospiraceae", "uncultured",
                     "Ruminococcus sp. ID1", "Ruminococcus")
A.slv_flt_GS_08 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_08[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_08[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_08[3])
A.slv_fix_GS_08 <- dplyr::mutate(
  A.slv_flt_GS_08,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_08[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_08[2],
                          values = A.slv_vec_GS_08[4]))

A.slv_vec_GS_09 <- c("Ruminococcaceae", "uncultured",
                     "Ruminococcaceae bacterium Marseille-P3738")
A.slv_flt_GS_09 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_09[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_09[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_09[3])
A.slv_fix_GS_09 <- dplyr::mutate(
  A.slv_flt_GS_09,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_09[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_09[3],
                            values = A.slv_vec_GS_09[2]))

A.slv_vec_GS_10 <- c("Lachnospiraceae", "uncultured",
                     "butyrate-producing bacterium L2-10")
A.slv_flt_GS_10 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_10[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_10[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_10[3])
A.slv_fix_GS_10 <- dplyr::mutate(
  A.slv_flt_GS_10,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_10[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_10[3],
                            values = A.slv_vec_GS_10[2]))

A.slv_vec_GS_11 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-41", "Clostridium")
A.slv_flt_GS_11 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_11[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_11[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_11[3])
A.slv_fix_GS_11 <- dplyr::mutate(
  A.slv_flt_GS_11,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_11[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_11[2],
                          values = A.slv_vec_GS_11[4]))

A.slv_vec_GS_12 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-54", "Clostridium")
A.slv_flt_GS_12 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_12[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_12[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_12[3])
A.slv_fix_GS_12 <- dplyr::mutate(
  A.slv_flt_GS_12,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_12[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_12[2],
                          values = A.slv_vec_GS_12[4]))

A.slv_vec_GS_13 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Clone-49", "Clostridium")
A.slv_flt_GS_13 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_13[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_13[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_13[3])
A.slv_fix_GS_13 <- dplyr::mutate(
  A.slv_flt_GS_13,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_13[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_13[2],
                          values = A.slv_vec_GS_13[4]))

A.slv_vec_GS_14 <- c("Lachnospiraceae", "uncultured",
                     "Clostridium sp. Culture-27", "Clostridium")
A.slv_flt_GS_14 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_14[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_14[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_14[3])
A.slv_fix_GS_14 <- dplyr::mutate(
  A.slv_flt_GS_14,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_14[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_14[2],
                          values = A.slv_vec_GS_14[4]))

A.slv_vec_GS_15 <- c("Barnesiellaceae", "uncultured",
                     "Bacteroides sp. Tilapia9", "Bacteroides")
A.slv_flt_GS_15 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_15[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_15[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_15[3])
A.slv_fix_GS_15 <- dplyr::mutate(
  A.slv_flt_GS_15,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_15[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_15[2],
                          values = A.slv_vec_GS_15[4]))

A.slv_vec_GS_16 <- c("Sphingobacteriaceae", "uncultured",
                     "Sphingobacterium jejuense", "Sphingobacterium")
A.slv_flt_GS_16 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_16[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_16[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_16[3])
A.slv_fix_GS_16 <- dplyr::mutate(
  A.slv_flt_GS_16,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_GS_16[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_GS_16[2],
                          values = A.slv_vec_GS_16[4]))

A.slv_vec_GS_17 <- c("Burkholderiaceae", "uncultured",
                     "Burkholderiales bacterium X4")
A.slv_flt_GS_17 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_17[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_17[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_17[3])
A.slv_fix_GS_17 <- dplyr::mutate(
  A.slv_flt_GS_17,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_17[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_17[3],
                            values = A.slv_vec_GS_17[2]))

A.slv_vec_GS_18 <- c("Burkholderiaceae", "uncultured",
                     "beta proteobacterium WX163")
A.slv_flt_GS_18 <- dplyr::filter(A.test_slv_trunc,
                                 .data[[A.slv_col[5]]] == A.slv_vec_GS_18[1],
                                 .data[[A.slv_col[6]]] == A.slv_vec_GS_18[2],
                                 .data[[A.slv_col[7]]] == A.slv_vec_GS_18[3])
A.slv_fix_GS_18 <- dplyr::mutate(
  A.slv_flt_GS_18,
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list = .data[[A.slv_col[5]]] == A.slv_vec_GS_18[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_GS_18[3],
                            values = A.slv_vec_GS_18[2]))

# Order-Family-Genus
A.slv_vec_OFG_01 <- c("Clostridiales", "Family XIII",
                      "Family XIII AD3011 group")
A.slv_flt_OFG_01 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_01[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_01[3])
A.slv_fix_OFG_01 <- dplyr::mutate(
  A.slv_flt_OFG_01,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_01[2],
                           values = paste(A.slv_vec_OFG_01[1],
                                          A.slv_vec_OFG_01[2],
                                          sep = " ")),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_01[3],
                          values = paste(A.slv_vec_OFG_01[1],
                                         A.slv_vec_OFG_01[3],
                                         sep = " ")))

A.slv_vec_OFG_02 <- c("Clostridiales", "Family XIII",
                      "Family XIII UCG-001")
A.slv_flt_OFG_02 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_02[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_02[3])
A.slv_fix_OFG_02 <- dplyr::mutate(
  A.slv_flt_OFG_02,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_02[2],
                           values = paste(A.slv_vec_OFG_02[1],
                                          A.slv_vec_OFG_02[2],
                                          sep = " ")),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_02[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_02[3],
                          values = paste(A.slv_vec_OFG_02[1],
                                         A.slv_vec_OFG_02[3],
                                         sep = " ")))

A.slv_vec_OFG_03 <- c("Rhodospirillales", "uncultured",
                      "Azospirillum sp. 47_25", "Rhodospirillaceae",
                      "Azospirillum")
A.slv_flt_OFG_03 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1],
                                  .data[[A.slv_col[5]]] == A.slv_vec_OFG_03[2],
                                  .data[[A.slv_col[6]]] == A.slv_vec_OFG_03[3])
A.slv_fix_OFG_03 <- dplyr::mutate(
  A.slv_flt_OFG_03,
  int.slv.Family = replace(x = .data[[A.slv_col[5]]],
                           list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1] &
                             .data[[A.slv_col[5]]] == A.slv_vec_OFG_03[2],
                           values = A.slv_vec_OFG_03[4]),
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[4]]] == A.slv_vec_OFG_03[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_OFG_03[3],
                          values = A.slv_vec_OFG_03[5]))

# Family-Genus-Species
A.slv_vec_FGS_01 <- c("Lachnospiraceae", "uncultured",
                      "intestinal bacterium CG19-1", "CG19-1")
A.slv_flt_FGS_01 <- dplyr::filter(A.test_slv_trunc,
                                  .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1],
                                  .data[[A.slv_col[6]]] == A.slv_vec_FGS_01[2],
                                  .data[[A.slv_col[7]]] == A.slv_vec_FGS_01[3])
A.slv_fix_FGS_01 <- dplyr::mutate(
  A.slv_flt_FGS_01,
  int.slv.Genus = replace(x = .data[[A.slv_col[6]]],
                          list = .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1] &
                            .data[[A.slv_col[6]]] == A.slv_vec_FGS_01[2],
                          values = paste(A.slv_vec_FGS_01[1],
                                         A.slv_vec_FGS_01[4],
                                         sep = " ")),
  int.slv.Species = replace(x = .data[[A.slv_col[7]]],
                            list =
                              .data[[A.slv_col[5]]] == A.slv_vec_FGS_01[1] &
                              .data[[A.slv_col[7]]] == A.slv_vec_FGS_01[3],
                            values = A.slv_vec_FGS_01[2]))

# combine the _flt dfs within groups and also combine the _fix dfs within groups
A.slv_flt_PC_all <- rbind(A.slv_flt_PC_01)
A.slv_fix_PC_all <- rbind(A.slv_fix_PC_01)

A.slv_flt_CO_all <- rbind(A.slv_flt_CO_01, A.slv_flt_CO_02, A.slv_flt_CO_03)
A.slv_fix_CO_all <- rbind(A.slv_fix_CO_01, A.slv_fix_CO_02, A.slv_fix_CO_03)

A.slv_flt_OF_all <- rbind(A.slv_flt_OF_01, A.slv_flt_OF_02, A.slv_flt_OF_03,
                          A.slv_flt_OF_04, A.slv_flt_OF_05, A.slv_flt_OF_06,
                          A.slv_flt_OF_07, A.slv_flt_OF_08, A.slv_flt_OF_09,
                          A.slv_flt_OF_10, A.slv_flt_OF_11, A.slv_flt_OF_12)
A.slv_fix_OF_all <- rbind(A.slv_fix_OF_01, A.slv_fix_OF_02, A.slv_fix_OF_03,
                          A.slv_fix_OF_04, A.slv_fix_OF_05, A.slv_fix_OF_06,
                          A.slv_fix_OF_07, A.slv_fix_OF_08, A.slv_fix_OF_09,
                          A.slv_fix_OF_10, A.slv_fix_OF_11, A.slv_fix_OF_12)

A.slv_flt_FG_all <- rbind(A.slv_flt_FG_01, A.slv_flt_FG_02, A.slv_flt_FG_03,
                          A.slv_flt_FG_04, A.slv_flt_FG_05, A.slv_flt_FG_06,
                          A.slv_flt_FG_07, A.slv_flt_FG_08, A.slv_flt_FG_09,
                          A.slv_flt_FG_10, A.slv_flt_FG_11, A.slv_flt_FG_12,
                          A.slv_flt_FG_13, A.slv_flt_FG_14, A.slv_flt_FG_15,
                          A.slv_flt_FG_16)
A.slv_fix_FG_all <- rbind(A.slv_fix_FG_01, A.slv_fix_FG_02, A.slv_fix_FG_03,
                          A.slv_fix_FG_04, A.slv_fix_FG_05, A.slv_fix_FG_06,
                          A.slv_fix_FG_07, A.slv_fix_FG_08, A.slv_fix_FG_09,
                          A.slv_fix_FG_10, A.slv_fix_FG_11, A.slv_fix_FG_12,
                          A.slv_fix_FG_13, A.slv_fix_FG_14, A.slv_fix_FG_15,
                          A.slv_fix_FG_16)

A.slv_flt_GS_all <- rbind(A.slv_flt_GS_01, A.slv_flt_GS_02, A.slv_flt_GS_03,
                          A.slv_flt_GS_04, A.slv_flt_GS_05, A.slv_flt_GS_06,
                          A.slv_flt_GS_07, A.slv_flt_GS_08, A.slv_flt_GS_09,
                          A.slv_flt_GS_10, A.slv_flt_GS_11, A.slv_flt_GS_12,
                          A.slv_flt_GS_13, A.slv_flt_GS_14, A.slv_flt_GS_15,
                          A.slv_flt_GS_16, A.slv_flt_GS_17, A.slv_flt_GS_18)
A.slv_fix_GS_all <- rbind(A.slv_fix_GS_01, A.slv_fix_GS_02, A.slv_fix_GS_03,
                          A.slv_fix_GS_04, A.slv_fix_GS_05, A.slv_fix_GS_06,
                          A.slv_fix_GS_07, A.slv_fix_GS_08, A.slv_fix_GS_09,
                          A.slv_fix_GS_10, A.slv_fix_GS_11, A.slv_fix_GS_12,
                          A.slv_fix_GS_13, A.slv_fix_GS_14, A.slv_fix_GS_15,
                          A.slv_fix_GS_16, A.slv_fix_GS_17, A.slv_fix_GS_18)

A.slv_flt_OFG_all <- rbind(A.slv_flt_OFG_01, A.slv_flt_OFG_02, A.slv_flt_OFG_03)
A.slv_fix_OFG_all <- rbind(A.slv_fix_OFG_01, A.slv_fix_OFG_02, A.slv_fix_OFG_03)

A.slv_flt_FGS_all <- rbind(A.slv_flt_FGS_01)
A.slv_fix_FGS_all <- rbind(A.slv_fix_FGS_01)

# then combine all groups together
A.slv_flt_all <- rbind(A.slv_flt_PC_all, A.slv_flt_CO_all, A.slv_flt_OF_all,
                       A.slv_flt_FG_all, A.slv_flt_GS_all, A.slv_flt_OFG_all,
                       A.slv_flt_FGS_all)

A.slv_fix_all <- rbind(A.slv_fix_PC_all, A.slv_fix_CO_all, A.slv_fix_OF_all,
                       A.slv_fix_FG_all, A.slv_fix_GS_all, A.slv_fix_OFG_all,
                       A.slv_fix_FGS_all)

# filter to remove the unfixed lineages and then add in the fixed lineages
A.int_slv_trunc_flt <- dplyr::anti_join(x = A.test_slv_trunc, y = A.slv_fix_all,
                                        by = "FeatureID")

A.int_slv_trunc_fix <- rbind(A.int_slv_trunc_flt, A.slv_fix_all)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost after "fixing" truncated lineages
A.test_slv_trunc_fix <- testif_fts_lost(raw_df = A.int_slv,
                                        test_df = A.int_slv_trunc_fix,
                                        test_type = "trunc_fix",
                                        tax_type = "SILVA",
                                        col_join = "FeatureID")
A.lost_slv_trunc_fix <- gather_lost_fts(raw_df = A.int_slv,
                                        test_df = A.int_slv_trunc_fix,
                                        test_type = "trunc",
                                        tax_type = "SILVA",
                                        col_join = "FeatureID")

# check to make sure there were no NAs produced in the fixed truncated df
# ** note for KDP: possible addition to future version of testif_fts_lost() ** #
A.int_slv_trunc_test_NAs <- dplyr::filter_all(A.int_slv_trunc_fix,
                                              dplyr::any_vars(is.na(.)))

LOGICAL_A.int_slv_trunc_test_NAs <- ifelse(nrow(A.int_slv_trunc_test_NAs) == 0,
                                           yes = T, no = F)

if (!isTRUE(LOGICAL_A.int_slv_trunc_test_NAs)) {
  stop("NA values produced during fix for SILVA truncated lineages;
       run: `View(A.int_slv_trunc_test_NAs)` to see rows with NA values")
}

### ************************************
### A - STEP 1c - format taxonomy tables: determine lowest assignment level ----
### ************************************

# ** note for KDP: AS determine lowest assignment level version 0.6 ** #

# create vectors useful to determine lowest assigned level for each feature
A.int_ggs_str <- "int.ggs."
A.int_slv_str <- "int.slv."
A.val_una <- "Unassigned"
A.flt_col <- "FeatureID"

# for SILVA assignments, create vectors of unique values...
# ... using rows in all of the truncated level cols
# NOTE: these are also useful when looking for assignments that need "fixed"
A.int_slv_K_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[1]])
A.int_slv_P_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[2]])
A.int_slv_C_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[3]])
A.int_slv_O_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[4]])
A.int_slv_F_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[5]])
A.int_slv_G_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[6]])
A.int_slv_S_unq <- unique(A.test_slv_trunc_fix[, A.slv_col[7]])

# create a vector of values we do not want in our lowest assignment column
# NOTE: uses partial string matching
A.slv_pattern <- paste(c("unass", "uncul", "unide", "metagen", "enrich"), 
                       collapse = "|")

# grep through the _unq dfs to create new vectors with the unwanted values
A.int_slv_K_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_K_unq,
                            ignore.case = T, value = T)
A.int_slv_K_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_K_unq,
                            ignore.case = T, value = T)
A.int_slv_P_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_P_unq,
                            ignore.case = T, value = T)
A.int_slv_C_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_C_unq,
                            ignore.case = T, value = T)
A.int_slv_O_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_O_unq,
                            ignore.case = T, value = T)
A.int_slv_F_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_F_unq,
                            ignore.case = T, value = T)
A.int_slv_G_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_G_unq,
                            ignore.case = T, value = T)
A.int_slv_S_unq_unw <- grep(pattern = A.slv_pattern, x = A.int_slv_S_unq,
                            ignore.case = T, value = T)

# combine the above vectors, pass unique() again, and collapse everything ...
# .... into a vector for use in grepl pattern matching below
A.val_slv <- paste(unique(c(A.int_slv_K_unq_unw, A.int_slv_P_unq_unw,
                            A.int_slv_C_unq_unw, A.int_slv_O_unq_unw,
                            A.int_slv_F_unq_unw, A.int_slv_G_unq_unw,
                            A.int_slv_S_unq_unw)), collapse = "|")

# now that we have our vector, process each of the truncated level cols...
# ... to find the lowest assignment level
# NOTE: any hashed lines encountered below had 0 obs. for these data... 
# ... i.e. there was no "Unassigned" for that respective level

# Unassigned Kingdom:
A.int_ggs_U_1 <- dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col),
                                  dplyr::all_vars(. == A.val_una))
A.int_ggs_U_2 <- dplyr::select(A.int_ggs_U_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[1]))
A.int_ggs_U_3 <- dplyr::rename(A.int_ggs_U_2, int.ggs.lws.txn = A.ggs_col[1])
A.int_ggs_U <- A.int_ggs_U_3
A.int_ggs_U$int.ggs.lws.lvl <- A.val_una

A.int_slv_U_1 <- dplyr::filter_at(A.test_slv_trunc_fix,
                                  dplyr::vars(A.slv_col),
                                  dplyr::all_vars(
                                    grepl(pattern = A.val_slv, x = .)))
A.int_slv_U_2 <- dplyr::select(A.int_slv_U_1,
                               dplyr::one_of(A.flt_col, A.slv_col[1]))
A.int_slv_U_3 <- dplyr::rename(A.int_slv_U_2, int.slv.lws.txn = A.slv_col[1])
A.int_slv_U <- A.int_slv_U_3
A.int_slv_U$int.slv.lws.lvl <- A.val_una

# Assigned Kingdom; Unassigned Phylum:
A.int_ggs_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[2:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_K_2 <- dplyr::select(A.int_ggs_K_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[1]))
A.int_ggs_K_3 <- dplyr::rename(A.int_ggs_K_2, int.ggs.lws.txn = A.ggs_col[1])
A.int_ggs_K <- A.int_ggs_K_3
A.int_ggs_K$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[1], A.int_ggs_str)[[1]][2]

A.int_slv_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[2:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_K_2 <- dplyr::select(A.int_slv_K_1,
                               dplyr::one_of(A.flt_col, A.slv_col[1]))
A.int_slv_K_3 <- dplyr::rename(A.int_slv_K_2, int.slv.lws.txn = A.slv_col[1])
A.int_slv_K <- A.int_slv_K_3
A.int_slv_K$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[1], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum; Unassigned Class:
A.int_ggs_P_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:2]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[3:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_P_2 <- dplyr::select(A.int_ggs_P_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[2]))
A.int_ggs_P_3 <- dplyr::rename(A.int_ggs_P_2, int.ggs.lws.txn = A.ggs_col[2])
A.int_ggs_P <- A.int_ggs_P_3
A.int_ggs_P$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[2], A.int_ggs_str)[[1]][2]

A.int_slv_P_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:2]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[3:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_P_2 <- dplyr::select(A.int_slv_P_1,
                               dplyr::one_of(A.flt_col, A.slv_col[2]))
A.int_slv_P_3 <- dplyr::rename(A.int_slv_P_2, int.slv.lws.txn = A.slv_col[2])
A.int_slv_P <- A.int_slv_P_3
A.int_slv_P$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[2], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class; Unassigned Order:
A.int_ggs_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:3]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[4:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_C_2 <- dplyr::select(A.int_ggs_C_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[3]))
A.int_ggs_C_3 <- dplyr::rename(A.int_ggs_C_2, int.ggs.lws.txn = A.ggs_col[3])
A.int_ggs_C <- A.int_ggs_C_3
A.int_ggs_C$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[3], A.int_ggs_str)[[1]][2]

A.int_slv_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:3]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[4:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_C_2 <- dplyr::select(A.int_slv_C_1,
                               dplyr::one_of(A.flt_col, A.slv_col[3]))
A.int_slv_C_3 <- dplyr::rename(A.int_slv_C_2, int.slv.lws.txn = A.slv_col[3])
A.int_slv_C <- A.int_slv_C_3
A.int_slv_C$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[3], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order; Unassigned Family:
A.int_ggs_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:4]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[5:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_O_2 <- dplyr::select(A.int_ggs_O_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[4]))
A.int_ggs_O_3 <- dplyr::rename(A.int_ggs_O_2, int.ggs.lws.txn = A.ggs_col[4])
A.int_ggs_O <- A.int_ggs_O_3
A.int_ggs_O$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[4], A.int_ggs_str)[[1]][2]

A.int_slv_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:4]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[5:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_O_2 <- dplyr::select(A.int_slv_O_1,
                               dplyr::one_of(A.flt_col, A.slv_col[4]))
A.int_slv_O_3 <- dplyr::rename(A.int_slv_O_2, int.slv.lws.txn = A.slv_col[4])
A.int_slv_O <- A.int_slv_O_3
A.int_slv_O$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[4], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family; Unassigned Genus:
A.int_ggs_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:5]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[6:7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_F_2 <- dplyr::select(A.int_ggs_F_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[5]))
A.int_ggs_F_3 <- dplyr::rename(A.int_ggs_F_2, int.ggs.lws.txn = A.ggs_col[5])
A.int_ggs_F <- A.int_ggs_F_3
A.int_ggs_F$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[5], A.int_ggs_str)[[1]][2]

A.int_slv_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:5]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[6:7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_F_2 <- dplyr::select(A.int_slv_F_1,
                               dplyr::one_of(A.flt_col, A.slv_col[5]))
A.int_slv_F_3 <- dplyr::rename(A.int_slv_F_2, int.slv.lws.txn = A.slv_col[5])
A.int_slv_F <- A.int_slv_F_3
A.int_slv_F$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[5], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family-Genus; Unassigned Species:
A.int_ggs_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[1:6]),
                   dplyr::all_vars(!. == A.val_una)),
  dplyr::filter_at(A.test_ggs_trunc, dplyr::vars(A.ggs_col[7]),
                   dplyr::all_vars(. == A.val_una)))
A.int_ggs_G_2 <- dplyr::select(A.int_ggs_G_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[6]))
A.int_ggs_G_3 <- dplyr::rename(A.int_ggs_G_2, int.ggs.lws.txn = A.ggs_col[6])
A.int_ggs_G <- A.int_ggs_G_3
A.int_ggs_G$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[6], A.int_ggs_str)[[1]][2]

A.int_slv_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[1:6]),
                   dplyr::all_vars(!grepl(pattern = A.val_slv, x = .))),
  dplyr::filter_at(A.test_slv_trunc_fix, dplyr::vars(A.slv_col[7]),
                   dplyr::all_vars(grepl(pattern = A.val_slv, x = .))))
A.int_slv_G_2 <- dplyr::select(A.int_slv_G_1,
                               dplyr::one_of(A.flt_col, A.slv_col[6]))
A.int_slv_G_3 <- dplyr::rename(A.int_slv_G_2, int.slv.lws.txn = A.slv_col[6])
A.int_slv_G <- A.int_slv_G_3
A.int_slv_G$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[6], A.int_slv_str)[[1]][2]

# Assigned Kingdom-Phylum-Class-Order-Family-Genus-Species:
A.int_ggs_S_1 <- dplyr::filter_at(A.test_ggs_trunc,
                                  dplyr::vars(A.ggs_col[1:7]),
                                  dplyr::all_vars(!. == A.val_una))
A.int_ggs_S_2 <- dplyr::select(A.int_ggs_S_1,
                               dplyr::one_of(A.flt_col, A.ggs_col[7]))
A.int_ggs_S_3 <- dplyr::rename(A.int_ggs_S_2, int.ggs.lws.txn = A.ggs_col[7])
A.int_ggs_S <- A.int_ggs_S_3
A.int_ggs_S$int.ggs.lws.lvl <- base::strsplit(
  A.ggs_col[7], A.int_ggs_str)[[1]][2]

A.int_slv_S_1 <- dplyr::filter_at(A.test_slv_trunc_fix,
                                  dplyr::vars(A.slv_col[1:7]),
                                  dplyr::all_vars(
                                    !grepl(pattern = A.val_slv, x = .)))
A.int_slv_S_2 <- dplyr::select(A.int_slv_S_1,
                               dplyr::one_of(A.flt_col, A.slv_col[7]))
A.int_slv_S_3 <- dplyr::rename(A.int_slv_S_2, int.slv.lws.txn = A.slv_col[7])
A.int_slv_S <- A.int_slv_S_3
A.int_slv_S$int.slv.lws.lvl <- base::strsplit(
  A.slv_col[7], A.int_slv_str)[[1]][2]

# combine the above dfs
A.int_ggs_lws <- rbind(
  A.int_ggs_U, 
  A.int_ggs_K, 
  A.int_ggs_P,
  A.int_ggs_C, 
  A.int_ggs_O, 
  A.int_ggs_F,
  A.int_ggs_G, 
  A.int_ggs_S)

A.int_slv_lws <- rbind(
  A.int_slv_U, 
  A.int_slv_K, 
  A.int_slv_P,
  A.int_slv_C, 
  A.int_slv_O, 
  A.int_slv_F,
  A.int_slv_G, 
  A.int_slv_S)

# run testif_fts_lost() and gather_lost_fts() functions to test for errors
# i.e. check that no features were lost in either of the lowest level dfs
A.test_ggs_lws <- testif_fts_lost(raw_df = A.int_ggs, 
                                  test_df = A.int_ggs_lws, 
                                  test_type = "lws", 
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.test_slv_lws <- testif_fts_lost(raw_df = A.int_slv, 
                                  test_df = A.int_slv_lws, 
                                  test_type = "lws",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

A.lost_ggs_lws <- gather_lost_fts(raw_df = A.int_ggs, 
                                  test_df = A.int_ggs_lws, 
                                  test_type = "lws", 
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.lost_slv_lws <- gather_lost_fts(raw_df = A.int_slv, 
                                  test_df = A.int_slv_lws, 
                                  test_type = "lws", 
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

# if no errors were produced, merge the lws dfs with the trunc or trunc_fix dfs
A.int_ggs_frm <- merge(A.test_ggs_trunc, A.test_ggs_lws, by = "FeatureID",
                       sort = F)
A.int_slv_frm <- merge(A.int_slv_trunc_fix, A.test_slv_lws, by = "FeatureID",
                       sort = F)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in the newly formatted dfs
A.test_ggs_frm <- testif_fts_lost(raw_df = A.int_ggs,
                                  test_df = A.int_ggs_frm,
                                  test_type = "frm",
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.test_slv_frm <- testif_fts_lost(raw_df = A.int_slv,
                                  test_df = A.int_slv_frm,
                                  test_type = "frm",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

A.lost_ggs_frm <- gather_lost_fts(raw_df = A.int_ggs,
                                  test_df = A.int_ggs_frm,
                                  test_type = "frm",
                                  tax_type = "Greengenes",
                                  col_join = "FeatureID")
A.lost_slv_frm <- gather_lost_fts(raw_df = A.int_slv,
                                  test_df = A.int_slv_frm,
                                  test_type = "frm",
                                  tax_type = "SILVA",
                                  col_join = "FeatureID")

# if no errors were produced, merge the formatted Greengenes and SILVA dfs
A.tax_mrg <- merge(A.test_ggs_frm, A.test_slv_frm, by = "FeatureID",
                   sort = F)

# run testif_fts_lost() and rather_lost_fts() functions to test for errors
# i.e. check that no features were lost in the merge
A.test_tax_mrg <- testif_fts_lost(raw_df = A.int_slv,
                                  test_df = A.tax_mrg,
                                  test_type = "mrg",
                                  tax_type = "merging",
                                  col_join = "FeatureID")

A.lost_tax_mrg <- gather_lost_fts(raw_df = A.int_slv,
                                  test_df = A.tax_mrg,
                                  test_type = "mrg",
                                  tax_type = "merging",
                                  col_join = "FeatureID")

### ************************************
### A - STEP  2 - create master table (EBTKS) ----
### ************************************

# ** note for KDP: AS create master table (EBTKS) version 0.7 ** #

# provide provenance for information gathering at end of STEP 2:
A.prov_secstep_AS2 <- "Section A - STEP 2"
A.prov_heading_AS2 <- "create master table (EBTKS)"
A.prov_output_obj_AS2 <- "A.EBTKS_abs_raw" # this object is output to the vault

## this step creates a master table containing:
# 'FeatureID' = MD5 sums generated in dada2 QIIME 2 feature table
# 'RepSeq' = representative DNA sequence for respective FeatureID
# 'int.ggs.tax' = taxonomic lineage with Greengenes naming convention
# 'int.ggs.cnf'= confidence in Greengenes taxonomic classification
# 'int.ggs.L' = taxon for respective level (L)
# 'int.ggs.lws.txn' = lowest assignment taxon by Greengenes
# 'int.ggs.lws.lvl' = lowest assignment level by Greengenes
# 'int.slv.tax' = taxonomic lineage with SILVA naming convention
# 'int.slv.cnf' = confidence in SILVA taxonomic classification
# 'int.slv.L' = taxon for respective level (L)
# 'int.slv.lws.txn' = lowest assignment taxon by SILVA
# 'int.slv.lws.lvl' = lowest assignment level by SILVA
# remaining columns = SampleID; row values = absolute counts of features

# read in dada2 feature table with absolute counts, and rename col #OTUID
A.raw_fts <- read.table(A.ifv_fts_tab, header = T, sep = "\t", as.is = T,
                        stringsAsFactors = F, check.names = F, skip = 1,
                        comment.char = "")

A.fts_abs <- dplyr::rename(A.raw_fts, FeatureID = "#OTU ID")

# read in the representative sequences fasta file as a data.frame (df)
A.rep_seq <- fasta_to_df(fasta_file = A.ifv_rep_seq, hdr_colname = "FeatureID",
                         seq_colname = "RepSeq")

# convert absolute feature counts to relative abundances
# NOTE: indexing with the numeric output of the which() function preserves... 
# ... column FeatureID no matter where it is located in the df
A.fts_rel <- A.fts_abs
A.fts_num <- which(names(A.fts_rel) == "FeatureID")
A.fts_rel[, -A.fts_num] <- lapply(A.fts_rel[, -A.fts_num],
                                  function(x) {x/sum(x)})

# optional, check conversion by running: colSums(A.fts_rel[, -A.fts_num])

# merge taxa, rep seqs, and fts dfs, & ensure that col FeatureID is a character
A.seq_tax <- merge(A.rep_seq, A.test_tax_mrg, by = "FeatureID", sort = F)
A.EBTKS_abs_raw <- merge(A.seq_tax, A.fts_abs, by = "FeatureID", sort = F)
A.EBTKS_rel_raw <- merge(A.seq_tax, A.fts_rel, by = "FeatureID", sort = F)

A.EBTKS_abs_raw$FeatureID <- as.character(A.EBTKS_abs_raw$FeatureID)
A.EBTKS_rel_raw$FeatureID <- as.character(A.EBTKS_rel_raw$FeatureID)

# ** note for KDP: something similar testif_fts_lost() might be needed here ** #

# gather info and provide some level of provenance outside of the script file
A.info_EBTKS_raw_AS2_a <- data.frame(
  "info" = "Number of samples",
  value = ncol(dplyr::select(A.EBTKS_abs_raw, -dplyr::one_of(com_col))),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)
A.info_EBTKS_raw_AS2_b <- data.frame(
  "info" = "Number of features",
  value = nrow(A.EBTKS_abs_raw),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)
A.info_EBTKS_raw_AS2_c <- data.frame(
  "info" = "Total frequency",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_raw, -dplyr::one_of(com_col)))),
  "object" = A.prov_output_obj_AS2, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS2, "heading" = A.prov_heading_AS2, 
  stringsAsFactors = F)

A.info_AS2 <- rbind(A.info_EBTKS_raw_AS2_a, A.info_EBTKS_raw_AS2_b,
                    A.info_EBTKS_raw_AS2_c)

### ************************************
### A - STEP 3a - process EBTKS: remove 'contaminants' ----
### ************************************

# ** note for KDP: AS remove 'contaminants' version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3a <- "Section A - STEP 3a"
A.prov_heading_AS3a <- "process EBTKS: remove 'contaminants'"
A.prov_object1_AS3a <- "A.ctm_abs"
A.prov_object2_AS3a <- "A.EBTKS_abs_pro_a"

# this sub-step uses taxonomic assignments from A.tax_mrg:
# to identify and remove features likely to be contamination ...
# ... e.g. chloroplast, mitohondria, eukaryotic, Unassigned Kingdom
A.clp_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "chloropl",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.clp_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "chloropl",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))
A.mtc_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "mitochon",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.mtc_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "mitochon",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))
A.euk_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "eukaryot",
                                                x = A.tax_mrg$int.ggs.tax,
                                                ignore.case = T))
A.euk_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "eukaryot",
                                                x = A.tax_mrg$int.slv.tax,
                                                ignore.case = T))

# careful here - CANNOT be partial match
A.unk_int_ggs <- dplyr::filter(A.tax_mrg, grepl(pattern = "Unassigned",
                                                x = A.tax_mrg$int.ggs.tax,
                                                fixed = T))
A.unk_int_slv <- dplyr::filter(A.tax_mrg, grepl(pattern = "Unassigned",
                                                x = A.tax_mrg$int.slv.tax,
                                                fixed = T))

# successively merge Greengenes and SILVA contaminant dfs
# unequal df lengths are handled by not specifying a 'by' variable
A.clp_int_mrg <- merge(A.clp_int_ggs, A.clp_int_slv, all = T, sort = F)
A.mtc_int_mrg <- merge(A.mtc_int_ggs, A.mtc_int_slv, all = T, sort = F)
A.euk_int_mrg <- merge(A.euk_int_ggs, A.euk_int_slv, all = T, sort = F)
A.unk_int_mrg <- merge(A.unk_int_ggs, A.unk_int_slv, all = T, sort = F)

# combine contaminant dfs
A.ctm <- rbind(A.clp_int_mrg, A.mtc_int_mrg, A.euk_int_mrg, A.unk_int_mrg)

# create ctm dfs with absolute or relative abundances
A.ctm_abs <- dplyr::semi_join(x = A.EBTKS_abs_raw, y = A.ctm, by = "FeatureID")
A.ctm_rel <- dplyr::semi_join(x = A.EBTKS_rel_raw, y = A.ctm, by = "FeatureID")

# remove contaminants from A.EBTKS_abs_raw
A.EBTKS_abs_pro_a <- dplyr::anti_join(x = A.EBTKS_abs_raw, y = A.ctm,
                                      by = "FeatureID")

### ************************************
### A - STEP 3b - process EBTKS: remove low feature count samples ----
### ************************************

# ** note for KDP: AS remove low feature count samples version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3b <- "Section A - STEP 3b"
A.prov_heading_AS3b <- "process EBTKS: remove low feature count samples"
A.prov_object1_AS3b <- "A.low_thrsh"
A.prov_object2_AS3b <- "A.low_fts_sID"
A.prov_object3_AS3b <- "A.EBTKS_abs_pro_b"

# this sub-step identifies and removes samples with:
# low total feature count for features not considered 'contaminants'
# NOTE: this step DOES NOT remove any protocol control samples... 
# ... these are expected to have low total feature counts

# after downstream processing, sample neg_cnt_12 was identified as causing ...
# ... an issue after clr transformation (i.e. returned negative values)
# this sample will be added to a separate object and removed in this step
A.negative_sID <- "neg_cnt_12"

# this process occurs as follows:
# (1) define the low total feature count threshold 
# (2) remove unneeded columns from A.EBTKS_abs_pro_a
# (3) sum absolute feature counts for each sample
# (4) create a vector of SampleIDs with low total feature count (to be removed)
# (5) ensure that SampleID is a character and not a factor
# (6) remove protocol controls from the sID vector
# (7) rename the sID vector
# (8) remove the identified samples from A.EBTKS_abs_pro_a

A.low_thrsh <- 400
A.EBTKS_abs_pro_frm <- dplyr::select(A.EBTKS_abs_pro_a, -dplyr::one_of(com_col))
A.EBTKS_abs_pro_fts_sum <- data.frame("SampleID" = names(A.EBTKS_abs_pro_frm),
                                      "Total" =  colSums(A.EBTKS_abs_pro_frm),
                                      row.names = NULL, stringsAsFactors = F)
A.low_fts_smp <- dplyr::filter(A.EBTKS_abs_pro_fts_sum, Total < A.low_thrsh)
A.low_fts_sID_0 <- as.character(A.low_fts_smp$SampleID)
A.low_fts_sID_1 <- A.low_fts_sID_0[!grepl(paste0("neg", collapse = "|"),
                                          A.low_fts_sID_0)]
A.low_fts_sID_2 <- A.low_fts_sID_1[!grepl(paste0("Swab", collapse = "|"),
                                          A.low_fts_sID_1)]
A.low_fts_sID_3 <- A.low_fts_sID_2[!grepl(paste0("Blank", collapse = "|"),
                                          A.low_fts_sID_2)]
A.low_fts_sID_4 <- A.low_fts_sID_3[!grepl(paste0("H2O", collapse = "|"),
                                          A.low_fts_sID_3)]
A.low_fts_sID_5 <- A.low_fts_sID_3[!grepl(paste0("cnt", collapse = "|"),
                                          A.low_fts_sID_3)]
A.low_fts_sID <- A.low_fts_sID_4
A.EBTKS_abs_pro_b <- dplyr::select(A.EBTKS_abs_pro_a, 
                                   -dplyr::one_of(A.negative_sID),
                                   -dplyr::one_of(A.low_fts_sID))

### ************************************
### A - STEP 3c - process EBTKS: remove high 'contaminant' count samples ----
### ************************************

# ** note for KDP: AS remove high 'contaminant' count samples version 0.5 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3c <- "Section A - STEP 3c"
A.prov_heading_AS3c <- "process EBTKS: remove high 'contaminant' count samples"
A.prov_object1_AS3c <- "A.hgh_thrsh"
A.prov_object2_AS3c <- "A.hgh_ctm_sID"

# this sub-step identifies and removes samples with:
# high total relative abundance for features considered 'contaminants'
# NOTE: this can remove protocol control samples

# this process occurs as follows:
# (1) define high total relative abundance threshold 
# (2) remove unneeded columns from A.ctm_rel
# (3) sum relative abundances for each sample
# (4) create a vector of SampleIDs with high contam. rel. abund. (to be removed)
# (5) ensure that SampleID is a character and not a factor
# (6) remove the identified samples from A.EBTKS_abs_pro_b

A.hgh_thrsh <- 0.05
A.ctm_rel_frm <- dplyr::select(A.ctm_rel, -dplyr::one_of(com_col))
A.ctm_rel_frm_sum <- data.frame("SampleID" = names(A.ctm_rel_frm),
                                "Total" = colSums(A.ctm_rel_frm),
                                row.names = NULL, stringsAsFactors = F)
A.hgh_ctm_smp <- dplyr::filter(A.ctm_rel_frm_sum, Total > A.hgh_thrsh)
A.hgh_ctm_sID <- as.character(A.hgh_ctm_smp$SampleID)
A.EBTKS_abs_pro_c <- dplyr::select(A.EBTKS_abs_pro_b,
                                   -dplyr::one_of(A.hgh_ctm_sID))

### ************************************
### A - STEP 3d - process EBTKS: remove zero count features ----
### ************************************

# ** note for KDP: AS remove zero count features version 0.7 ** #

# provide provenance for information gathering in STEP 3i:
A.prov_secstep_AS3d <- "Section A - STEP 3d"
A.prov_heading_AS3d <- "process EBTKS: remove zero count features"
A.prov_output_obj_AS3d <- "A.EBTKS_abs_pro" # this object is output to the vault

# this sub-step identifies and removes features with:
# a total count of less than 2* across all samples^
# *threshold is in accordance with the dada2 feature table produced in QIIME 2
# ^this can & often does occur when samples are removed from a feature table ...
# ... which is a possibility after running through A - STEP 3c above
# zero count features should be removed prior to any downstream analysis

# this process occurs as follows:
# (1) create a copy of A.EBTKS_abs_pro_c to avoid overwriting data
# (2) convert column FeatureID into row.names
# (3) remove unneeded columns
# (4) sum absolute counts for each feature across all samples
# (5) create a data.frame of features that will be removed
# (6) ensure that FeatureID is a character and not a factor
# (7) anti join A.EBTKS_abs_pro_c to remove all features present in the rmv df
# NOTE: the naming convention where the sub-step's letter is appended ...
# ... to the end of the processed data.frame is dropped in this sub-step

A.EBTKS_abs_pro_1 <- A.EBTKS_abs_pro_c
row.names(A.EBTKS_abs_pro_1) <- A.EBTKS_abs_pro_1$FeatureID
A.EBTKS_abs_pro_2 <- dplyr::select(A.EBTKS_abs_pro_1, -dplyr::one_of(com_col))
A.EBTKS_abs_pro_sum <- data.frame("FeatureID" = row.names(A.EBTKS_abs_pro_2),
                                  "FeatureTotal" = rowSums(A.EBTKS_abs_pro_2),
                                  row.names = NULL, stringsAsFactors = F)
A.EBTKS_abs_pro_rmv <- dplyr::filter(A.EBTKS_abs_pro_sum, FeatureTotal < 2)
A.EBTKS_abs_pro_rmv$FeatureID <- as.character(A.EBTKS_abs_pro_rmv$FeatureID)
A.EBTKS_abs_pro <- dplyr::anti_join(x = A.EBTKS_abs_pro_c, 
                                    y = A.EBTKS_abs_pro_rmv, by = "FeatureID")

### ************************************
### A - STEP 3i - information gathering and provenance ----
### ************************************

# ** note for KDP: AS information gathering and provenance version 0.2 ** #

## A - STEP 3a:
A.info_ctm_abs_AS3a_a <- data.frame(
  "info" = "Number of features considered 'contaminants'",
  value = nrow(A.ctm_abs),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_b <- data.frame(
  "info" = "Total frequency of 'contaminants'",
  value = sum(colSums(dplyr::select(A.ctm_abs, -dplyr::one_of(com_col)))),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_c <- data.frame(
  "info" = "Number of samples 'contaminants' observed in",
  value = length(which(colSums(
    dplyr::select(A.ctm_abs, -dplyr::one_of(com_col))) > 0)),
  "object" = A.prov_object1_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_ctm_abs_AS3a_d_0 <- data.frame(
  value = which(colSums(dplyr::select(A.ctm_abs, -dplyr::one_of(com_col))) > 0))
A.info_ctm_abs_AS3a_d_1 <- A.info_ctm_abs_AS3a_d_0
A.info_ctm_abs_AS3a_d_1$info <- paste(row.names(A.info_ctm_abs_AS3a_d_1), 
                                      ": Total frequency of 'contaminants'", 
                                      sep = "")
A.info_ctm_abs_AS3a_d_2 <- A.info_ctm_abs_AS3a_d_1
A.info_ctm_abs_AS3a_d_2$object <- A.prov_object1_AS3a
A.info_ctm_abs_AS3a_d_2$script <- name_scrpt
A.info_ctm_abs_AS3a_d_2$section <- A.prov_secstep_AS3a
A.info_ctm_abs_AS3a_d_2$heading <- A.prov_heading_AS3a
A.info_ctm_abs_AS3a_d <- A.info_ctm_abs_AS3a_d_2
row.names(A.info_ctm_abs_AS3a_d) <- 1:nrow(A.info_ctm_abs_AS3a_d)

A.info_ctm_abs_AS3a <- rbind(A.info_ctm_abs_AS3a_a, A.info_ctm_abs_AS3a_b, 
                             A.info_ctm_abs_AS3a_c, A.info_ctm_abs_AS3a_d)

A.info_EBTKS_abs_pro_a_AS3a_a <- data.frame(
  "info" = "Number of samples after 'contaminant' removal",
  value = ncol(dplyr::select(A.EBTKS_abs_pro_a, -dplyr::one_of(com_col))),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_a_AS3a_b <- data.frame(
  "info" = "Number of features after 'contaminant' removal",
  value = nrow(A.EBTKS_abs_pro_a),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_a_AS3a_c <- data.frame(
  "info" = "Total frequency after 'contaminant' removal",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro_a, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_object2_AS3a, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3a, "heading" = A.prov_heading_AS3a, 
  stringsAsFactors = F)

A.info_EBTKS_abs_pro_a_AS3a <- rbind(A.info_EBTKS_abs_pro_a_AS3a_a, 
                                     A.info_EBTKS_abs_pro_a_AS3a_b,
                                     A.info_EBTKS_abs_pro_a_AS3a_c)

A.info_AS3a <- rbind(A.info_ctm_abs_AS3a, A.info_EBTKS_abs_pro_a_AS3a)

## A - STEP 3b:
A.info_low_fts_AS3b_a <- data.frame(
  "info" = "low total feature count threshold",
  value = paste("< ", A.low_thrsh, sep = ""),
  "object" = A.prov_object1_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_AS3b_b <- data.frame(
  "info" = "Number of samples below threshold",
  value = length(A.low_fts_sID),
  "object" = A.prov_object2_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_AS3b_c <- data.frame(
  "info" = "SampleID below threshold",
  value = A.low_fts_sID,
  "object" = A.prov_object2_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_low_fts_AS3b_d <- data.frame(
  "info" = "negative control SampleID removed",
  value = A.negative_sID,
  "object" = A.prov_object2_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)

A.info_low_fts_smp_AS3b <- rbind(A.info_low_fts_AS3b_a, A.info_low_fts_AS3b_b,
                                 A.info_low_fts_AS3b_c, A.info_low_fts_AS3b_d)

A.info_EBTKS_abs_pro_b_AS3b_a <- data.frame(
  "info" = "Number of samples after removal of low feature count samples",
  value = ncol(dplyr::select(A.EBTKS_abs_pro_b, -dplyr::one_of(com_col))),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_b_AS3b_b <- data.frame(
  "info" = "Number of features after removal of low feature count samples",
  value = nrow(A.EBTKS_abs_pro_b),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_b_AS3b_c <- data.frame(
  "info" = "Total frequency after removal of low feature count samples",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro_b, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_object3_AS3b, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3b, "heading" = A.prov_heading_AS3b, 
  stringsAsFactors = F)

A.info_EBTKS_abs_pro_b_AS3b <- rbind(A.info_EBTKS_abs_pro_b_AS3b_a, 
                                     A.info_EBTKS_abs_pro_b_AS3b_b,
                                     A.info_EBTKS_abs_pro_b_AS3b_c)

A.info_AS3b <- rbind(A.info_low_fts_smp_AS3b, A.info_EBTKS_abs_pro_b_AS3b)

## A - STEP 3c:
A.info_hgh_ctm_AS3c_a <- data.frame(
  "info" = "high 'contaminant' threshold",
  value = paste("> ", A.hgh_thrsh * 100, "%", sep = ""),
  "object" = A.prov_object1_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_hgh_ctm_AS3c_b <- data.frame(
  "info" = "Number of samples above threshold",
  value = length(A.hgh_ctm_sID),
  "object" = A.prov_object2_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_hgh_ctm_AS3c_c <- data.frame(
  "info" = "SampleID below threshold",
  value = A.hgh_ctm_sID,
  "object" = A.prov_object2_AS3c, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3c, "heading" = A.prov_heading_AS3c, 
  stringsAsFactors = F)
A.info_AS3c <- rbind(A.info_hgh_ctm_AS3c_a, A.info_hgh_ctm_AS3c_b,
                     A.info_hgh_ctm_AS3c_c)

# NOTE: info related to ...
# "Number of samples after removal of high 'contaminant' samples"
# "Number of features after removal of high 'contaminant' samples"
# "Total frequency after removal of high 'contaminant' samples"
# ... SHOULD NOT be calculated using the A.EBTKS_abs_pro_c df
# to understand why, read the comments at the the beginning of STEP 3d ...
# ... regarding what can happen when a sample is removed
# also, the info mentioned above is identical to the info gathered for STEP 3d

## A - STEP 3d:
A.info_EBTKS_abs_pro_AS3d_a <- data.frame(
  "info" = "Number of samples after 'process EBTKS' steps",
  value = ncol(dplyr::select(A.EBTKS_abs_pro, -dplyr::one_of(com_col))),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_AS3d_b <- data.frame(
  "info" = "Number of features after 'process EBTKS' steps",
  value = nrow(A.EBTKS_abs_pro),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)
A.info_EBTKS_abs_pro_AS3d_c <- data.frame(
  "info" = "Total frequency after 'process EBTKS' steps",
  value = sum(colSums(dplyr::select(A.EBTKS_abs_pro, 
                                    -dplyr::one_of(com_col)))),
  "object" = A.prov_output_obj_AS3d, "script" = name_scrpt, 
  "section" = A.prov_secstep_AS3d, "heading" = A.prov_heading_AS3d, 
  stringsAsFactors = F)

A.info_AS3d <- rbind(A.info_EBTKS_abs_pro_AS3d_a, A.info_EBTKS_abs_pro_AS3d_b,
                     A.info_EBTKS_abs_pro_AS3d_c)

# combine all A.info_ dfs created up to this point
A.info_AS3 <- rbind(A.info_AS3a, A.info_AS3b, A.info_AS3c, A.info_AS3d)

### ************************************
### A - WRITE OUTPUTS ----
### ************************************

# rbind A.info data.frames together
A.info <- rbind(A.info_AS2, A.info_AS3)

# provenance for outputs to the vault
A.prov_output1 <- data.frame("info" = "provenance for output",
                             "path" = A.ofv_EBTKS_raw,
                             "object" = A.prov_output_obj_AS2,
                             "script" = name_scrpt,
                             "section" = A.prov_secstep_AS2,
                             "heading" = A.prov_heading_AS2,
                             stringsAsFactors = F)
A.prov_output2 <- data.frame("info" = "provenance for output",
                             "path" = A.ofv_EBTKS_pro,
                             "object" = A.prov_output_obj_AS3d,
                             "script" = name_scrpt,
                             "section" = A.prov_secstep_AS3d,
                             "heading" = A.prov_heading_AS3d,
                             stringsAsFactors = F)
A.prov <- rbind(A.prov_output1, A.prov_output2)

# outputs to the vault
write.table(sep = "\t", row.names = F, 
            x = A.EBTKS_abs_raw, file = A.ofv_EBTKS_raw)
write.table(sep = "\t", row.names = F, 
            x = A.EBTKS_abs_pro, file = A.ofv_EBTKS_pro)

write.table(sep = "\t", row.names = F, x = A.info, file = A.ofv_info)
write.table(sep = "\t", row.names = F, x = A.prov, file = A.ofv_prov)

A.obj <- ls(pattern = "A.")
A.lst <- c(A.obj[grep(pattern = "A.", x = A.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS)
save(list = A.lst, file = A.ofv_wksp)

### ************************************
### S - STEP 1 - subset EBTKS: split by Cohort and Type ----
### ************************************

# ** note for KDP: SS subset EBTKS - version 0.2 ** #

# this section subsets A.EBTKS_abs_pro by Cohort and Type ...
# ... where cohort represents healthy or AOM_DSS treated mice... 
# ... and Type represents type of sample (i.e. location) cecum/colon/feces

# provide provenance for information gathering in STEP 2:
S.prov_secstep_SS1 <- "Section S - STEP 1"
S.prov_heading_SS1 <- "subset EBTKS: split by Cohort and Type"
S.prov_object_SecA <- "A.EBTKS_abs_pro"
S.prov_info <- paste("groups isolated from", S.prov_object_SecA, sep = " ")

# ** note for KDP: SS subset EBTKS - version 0.1 ** #

# NOTE: if the environment is empty; there are some requirements:
# Section S requires access to objects in PREFACE; Section A
# ** note for KDP:
# the choose your own adventure function that allows sections to standalone ...
# ... or stitches them together is still under construction ** #

## this section subsets A.EBTKS_abs_pro by Type and DietGroup
# abbreviations are as follows:
# _aom = AOM_DSS treated mice
# _hth = healthy mice
# _EPDf = cEcum; Proximal colon; Distal colon; feces
# _EPD = cEcum; Proximal colon; Distal colon
# _PD = Proximal colon; Distal colon
# _E = cEcum
# _P = Proximal colon
# _D = Distal colon

# create new versions of the objects needed from section A ...
# ... and create a vector naming those objects (used when saving the workspace)
S.EBTKS_abs_pro <- A.EBTKS_abs_pro
S.low_fts_sID <- A.low_fts_sID
S.negative_sID <- A.negative_sID
S.hgh_ctm_sID <- A.hgh_ctm_sID
S.obj_from_A <- c("A.EBTKS_abs_pro", "A.low_fts_sID", "A.negative_sID", 
                  "A.hgh_ctm_sID")

# the subsetting process occurs as follows:
# first: remove samples from smp_dat that were removed during EBTKS processing
# then: isolate samples; remove zero count features; track provenance
# rinse and repeat 'then:' for the other relevant subsets
# NOTE: removal of zero count feaures follows method described in A - STEP 3d

# first:

S.smp_dat_pro_b0 <- dplyr::filter(smp_dat, !SampleID %in% S.negative_sID)
S.smp_dat_pro_b <- dplyr::filter(S.smp_dat_pro_b0, !SampleID %in% S.low_fts_sID)
S.smp_dat_pro_c <- dplyr::filter(S.smp_dat_pro_b, !SampleID %in% S.hgh_ctm_sID)
S.smp_dat_pro <- S.smp_dat_pro_c

S.smp_dat_pro_aom_EPDf <- dplyr::filter(S.smp_dat_pro, 
                                        Cohort == "AOM_DSS", 
                                        Type == "cEcum_material" |
                                          Type == "colon_Proximal" |
                                          Type == "colon_Distal" |
                                          Type == "feces")
S.smp_dat_pro_aom_EPD <- dplyr::filter(S.smp_dat_pro, 
                                       Cohort == "AOM_DSS", 
                                       Type == "cEcum_material" |
                                         Type == "colon_Proximal" |
                                         Type == "colon_Distal")
S.smp_dat_pro_aom_PD <- dplyr::filter(S.smp_dat_pro, 
                                      Cohort == "AOM_DSS", 
                                      Type == "colon_Proximal" |
                                        Type == "colon_Distal")
S.smp_dat_pro_aom_E <- dplyr::filter(S.smp_dat_pro, 
                                     Cohort == "AOM_DSS", 
                                     Type == "cEcum_material")
S.smp_dat_pro_aom_P <- dplyr::filter(S.smp_dat_pro, 
                                     Cohort == "AOM_DSS", 
                                     Type == "colon_Proximal")
S.smp_dat_pro_aom_D <- dplyr::filter(S.smp_dat_pro, 
                                     Cohort == "AOM_DSS", 
                                     Type == "colon_Distal")

S.sID_aom_EPDf <- as.vector(S.smp_dat_pro_aom_EPDf$SampleID)
S.sID_aom_EPD <- as.vector(S.smp_dat_pro_aom_EPD$SampleID)
S.sID_aom_PD <- as.vector(S.smp_dat_pro_aom_PD$SampleID)
S.sID_aom_E <- as.vector(S.smp_dat_pro_aom_E$SampleID)
S.sID_aom_P <- as.vector(S.smp_dat_pro_aom_P$SampleID)
S.sID_aom_D <- as.vector(S.smp_dat_pro_aom_D$SampleID)

# then:

# aom_EPDf = AOM_DSS group cEcum; Proximal colon; Distal colon; feces
S.abs_aom_EPDf_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                                  dplyr::one_of(S.sID_aom_EPDf))
S.abs_aom_EPDf_1 <- S.abs_aom_EPDf_0
row.names(S.abs_aom_EPDf_1) <- S.abs_aom_EPDf_1$FeatureID
S.abs_aom_EPDf_2 <- dplyr::select(S.abs_aom_EPDf_1, -dplyr::one_of(com_col))
S.abs_aom_EPDf_sum <- data.frame("FeatureID" = row.names(S.abs_aom_EPDf_2),
                                 "FeatureTotal" = rowSums(S.abs_aom_EPDf_2),
                                 row.names = NULL)
S.abs_aom_EPDf_rmv <- dplyr::filter(S.abs_aom_EPDf_sum, FeatureTotal < 2)
S.abs_aom_EPDf_rmv$FeatureID <- as.character(S.abs_aom_EPDf_rmv$FeatureID)
S.abs_aom_EPDf <- dplyr::anti_join(x = S.abs_aom_EPDf_0, y = S.abs_aom_EPDf_rmv, 
                                   by = "FeatureID")
S.info_abs_aom_EPDf <- data.frame(
  "object" = "S.abs_aom_EPDf",
  "value" = "AOM_DSS group cEcum; Proximal colon; Distal colon; feces", 
  stringsAsFactors = F)

# aom_EPD = AOM_DSS group cEcum; Proximal colon; Distal colon
S.abs_aom_EPD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                                 dplyr::one_of(S.sID_aom_EPD))
S.abs_aom_EPD_1 <- S.abs_aom_EPD_0
row.names(S.abs_aom_EPD_1) <- S.abs_aom_EPD_1$FeatureID
S.abs_aom_EPD_2 <- dplyr::select(S.abs_aom_EPD_1, -dplyr::one_of(com_col))
S.abs_aom_EPD_sum <- data.frame("FeatureID" = row.names(S.abs_aom_EPD_2),
                                "FeatureTotal" = rowSums(S.abs_aom_EPD_2),
                                row.names = NULL)
S.abs_aom_EPD_rmv <- dplyr::filter(S.abs_aom_EPD_sum, FeatureTotal < 2)
S.abs_aom_EPD_rmv$FeatureID <- as.character(S.abs_aom_EPD_rmv$FeatureID)
S.abs_aom_EPD <- dplyr::anti_join(x = S.abs_aom_EPD_0, y = S.abs_aom_EPD_rmv, 
                                  by = "FeatureID")
S.info_abs_aom_EPD <- data.frame(
  "object" = "S.abs_aom_EPD",
  "value" = "AOM_DSS group cEcum; Proximal colon; Distal colon", 
  stringsAsFactors = F)

# aom_PD = AOM_DSS group Proximal colon; Distal colon
S.abs_aom_PD_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                                dplyr::one_of(S.sID_aom_PD))
S.abs_aom_PD_1 <- S.abs_aom_PD_0
row.names(S.abs_aom_PD_1) <- S.abs_aom_PD_1$FeatureID
S.abs_aom_PD_2 <- dplyr::select(S.abs_aom_PD_1, -dplyr::one_of(com_col))
S.abs_aom_PD_sum <- data.frame("FeatureID" = row.names(S.abs_aom_PD_2),
                               "FeatureTotal" = rowSums(S.abs_aom_PD_2),
                               row.names = NULL)
S.abs_aom_PD_rmv <- dplyr::filter(S.abs_aom_PD_sum, FeatureTotal < 2)
S.abs_aom_PD_rmv$FeatureID <- as.character(S.abs_aom_PD_rmv$FeatureID)
S.abs_aom_PD <- dplyr::anti_join(x = S.abs_aom_PD_0, y = S.abs_aom_PD_rmv, 
                                 by = "FeatureID")
S.info_abs_aom_PD <- data.frame(
  "object" = "S.abs_aom_PD",
  "value" = "AOM_DSS group Proximal colon; Distal colon", 
  stringsAsFactors = F)

# aom_E = AOM_DSS group cEcum
S.abs_aom_E_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                               dplyr::one_of(S.sID_aom_E))
S.abs_aom_E_1 <- S.abs_aom_E_0
row.names(S.abs_aom_E_1) <- S.abs_aom_E_1$FeatureID
S.abs_aom_E_2 <- dplyr::select(S.abs_aom_E_1, -dplyr::one_of(com_col))
S.abs_aom_E_sum <- data.frame("FeatureID" = row.names(S.abs_aom_E_2),
                              "FeatureTotal" = rowSums(S.abs_aom_E_2),
                              row.names = NULL)
S.abs_aom_E_rmv <- dplyr::filter(S.abs_aom_E_sum, FeatureTotal < 2)
S.abs_aom_E_rmv$FeatureID <- as.character(S.abs_aom_E_rmv$FeatureID)
S.abs_aom_E <- dplyr::anti_join(x = S.abs_aom_E_0, y = S.abs_aom_E_rmv, 
                                by = "FeatureID")
S.info_abs_aom_E <- data.frame(
  "object" = "S.abs_aom_E",
  "value" = "AOM_DSS group cEcum", 
  stringsAsFactors = F)

# aom_P = AOM_DSS group Proximal colon
S.abs_aom_P_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                               dplyr::one_of(S.sID_aom_P))
S.abs_aom_P_1 <- S.abs_aom_P_0
row.names(S.abs_aom_P_1) <- S.abs_aom_P_1$FeatureID
S.abs_aom_P_2 <- dplyr::select(S.abs_aom_P_1, -dplyr::one_of(com_col))
S.abs_aom_P_sum <- data.frame("FeatureID" = row.names(S.abs_aom_P_2),
                              "FeatureTotal" = rowSums(S.abs_aom_P_2),
                              row.names = NULL)
S.abs_aom_P_rmv <- dplyr::filter(S.abs_aom_P_sum, FeatureTotal < 2)
S.abs_aom_P_rmv$FeatureID <- as.character(S.abs_aom_P_rmv$FeatureID)
S.abs_aom_P <- dplyr::anti_join(x = S.abs_aom_P_0, y = S.abs_aom_P_rmv, 
                                by = "FeatureID")
S.info_abs_aom_P <- data.frame(
  "object" = "S.abs_aom_P",
  "value" = "AOM_DSS group Proximal colon", 
  stringsAsFactors = F)

# aom_D = AOM_DSS group Distal colon
S.abs_aom_D_0 <- dplyr::select(S.EBTKS_abs_pro, dplyr::one_of(com_col),
                               dplyr::one_of(S.sID_aom_D))
S.abs_aom_D_1 <- S.abs_aom_D_0
row.names(S.abs_aom_D_1) <- S.abs_aom_D_1$FeatureID
S.abs_aom_D_2 <- dplyr::select(S.abs_aom_D_1, -dplyr::one_of(com_col))
S.abs_aom_D_sum <- data.frame("FeatureID" = row.names(S.abs_aom_D_2),
                              "FeatureTotal" = rowSums(S.abs_aom_D_2),
                              row.names = NULL)
S.abs_aom_D_rmv <- dplyr::filter(S.abs_aom_D_sum, FeatureTotal < 2)
S.abs_aom_D_rmv$FeatureID <- as.character(S.abs_aom_D_rmv$FeatureID)
S.abs_aom_D <- dplyr::anti_join(x = S.abs_aom_D_0, y = S.abs_aom_D_rmv, 
                                by = "FeatureID")
S.info_abs_aom_D <- data.frame(
  "object" = "S.abs_aom_D",
  "value" = "AOM_DSS group Distal colon", 
  stringsAsFactors = F)

### ************************************
### S - STEP 2 - information gathering and provenance ----
### ************************************

# ** note for KDP: SS information gathering and provenance version 0.1 ** #

# (1) bind S.info data.frames together
# (2) add info/script/section/heading columns
# (3) reorder columns

S.info_0 <- rbind(S.info_abs_aom_EPDf,
                  S.info_abs_aom_EPD,
                  S.info_abs_aom_PD,
                  S.info_abs_aom_E,
                  S.info_abs_aom_P,
                  S.info_abs_aom_D)

S.info_1 <- S.info_0
S.info_1$info <- S.prov_info
S.info_1$script <- name_scrpt
S.info_1$section <- S.prov_secstep_SS1
S.info_1$heading <- S.prov_heading_SS1

S.info <- dplyr::select(S.info_1, info, value, object, script, section, heading)

### ************************************
### S - WRITE OUTPUTS ----
### ************************************

# outputs to the vault
write.table(sep = "\t", row.names = F, x = S.info, file = S.ofv_info)

S.obj <- ls(pattern = "S.")
S.lst <- c(S.obj[grep(pattern = "S.", x = S.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS, S.obj_from_A)
save(list = S.lst, file = S.ofv_wksp)

# ### ************************************
# ### T - STEP 1 - create a rooted phylogenetic tree from processed EBTKS ----
# ### ************************************
#
# # ** note for KDP: Section T version 0.6 ** #
# 
# # NOTE: tree construction follows the code structure in:
# # Microbiome Data Analysis: from raw reads to community analyses.
# # created by Benjamin J Callahan, Kris Sankaran, Julia A Fukuyama,
# # Paul Joey McMurdie, and Susan P Holmes
# 
# # NOTE: if the environment is empty; there are some requirements:
# # Section T requires access to objects in PREFACE; Section A
# 
# # provide provenance for information gathering at end of section:
# T.prov_secstep_TS1 <- "Section T - STEP 1"
# T.prov_heading_TS1 <- "create a rooted phylogenetic tree from processed EBTKS"
# T.prov_output_obj_TS1 <- "T.EBTKS_pro_tree" # this object is output to the vault
# 
# # create a new version of the object needed from section A ...
# # ... and create a vector naming that object (used when saving the workspace)
# T.EBTKS_abs_pro <- A.EBTKS_abs_pro
# T.obj_from_A <- "A.EBTKS_abs_pro"
# 
# # the tree construction process occurs as follows:
# # (1) create a vector of RepSeqs
# # (2) create a copy to avoid overwriting the original
# # (3) assign RepSeqs as explicit names
# # (4) align RepSeqs
# # (5) convert alignment to class 'phyDat'
# # (6) compute pairwise distances from DNA sequences
# # (7) construct a neighbor-joining tree
# # (8:10) fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree
# # (11) extract the unrooted tree
# # (12) create a midpoint rooted tree
# # (13) and then rename to include the EBTKS_pro naming convention
# # NOTE: GTR+G+I = Generalized time-reversible with Gamma rate variation
# # runtime (1:12): ~2 hours
# 
# T.RepSeq_0 <- T.EBTKS_abs_pro$RepSeq
# T.RepSeq_1 <- T.RepSeq_0
# names(T.RepSeq_1) <- T.RepSeq_1
# T.RepSeq_algn_0 <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(T.RepSeq_1),
#                                        anchor = NA, verbose = F,
#                                        processors = NULL)
# T.RepSeq_algn_1 <- phangorn::phyDat(as(T.RepSeq_algn_0, "matrix"), type = "DNA")
# T.dist_ml <- phangorn::dist.ml(T.RepSeq_algn_1, processors = NULL)
# T.tree_NJ <- phangorn::NJ(T.dist_ml)
# T.fit = phangorn::pml(T.tree_NJ, data = T.RepSeq_algn_1, processors = NULL)
# T.fitGTR <- update(T.fit, k = 4, inv = 0.2)
# T.fitGTR <- phangorn::optim.pml(T.fitGTR, model = "GTR", optInv = T,
#                                 optGamma = T, rearrangement = "stochastic",
#                                 control = phangorn::pml.control(trace = 0))
# T.tree_GTR <- T.fitGTR$tree
# T.tree_root_GTR <- phangorn::midpoint(tree = T.tree_GTR)
# T.EBTKS_pro_tree <- T.tree_root_GTR
# 
# ### ************************************
# ### T - WRITE OUTPUTS ----
# ### ************************************
# 
# # provenance for output to the vault
# T.prov <- data.frame("info" = "provenance for output",
#                      "path" = T.ofv_EBTKS_pro_tree,
#                      "object" = T.prov_output_obj_TS1,
#                      "script" = name_scrpt,
#                      "section" = T.prov_secstep_TS1,
#                      "heading" = T.prov_heading_TS1,
#                      stringsAsFactors = F)
# 
# # outputs to the vault
# ape::write.tree(phy = T.EBTKS_pro_tree, file = T.ofv_EBTKS_pro_tree)
# write.table(sep = "\t", row.names = F, x = T.prov, file = T.ofv_prov)
# 
# T.obj <- ls(pattern = "T.")
# T.lst <- c(T.obj[grep(pattern = "T.", x = T.obj, ignore.case = F, fixed = T)],
#            PREFACE.lst, T.obj_from_A)
# save(list = T.lst, file = T.ofv_wksp)

### ************************************
### E - FUNCTION ----
### ************************************

# create a vector naming all Section B functions (used when saving workspace)
E.function <- "format_wi_kw"

# ** note for KDP: version 0.1 ** #
# function takexs an input of one of the following classes:
# "htest" from kruskal.test()
# "pairwise.htest" from pairwise.wilcox.test()
# "PMCMR" from PMCMR::posthoc.kruskal.dunn.test()
# and returns a formatted results data.frame

format_wi_kw <- function(data, digits, cutpoints, symbols, names_P = c("","")) {
  
  # internal checks for correct input types
  if (!inherits(data, "htest") && !inherits(data, "pairwise.htest")
      && !inherits(data, "PMCMR")) {
    stop("data must be of class 'htest' or 'pairwise.htest' or 'PMCMR'")
  }
  if (!inherits(digits, "numeric")) {
    stop("input for digits must be class numeric")
  }
  if (!inherits(cutpoints, "numeric")) {
    stop("input for cutpoints must be class numeric")
  }
  if (!inherits(symbols, "character")) {
    stop("input for symbol must be class character")
  }
  if (!inherits(names_P, "character")) {
    stop("input for names_P must be class character")
  }
  if (!length(names_P) == 2) {
    stop("input for names_P must be vector of length 2")
  }
  
  # input data are from kruskal.test()
  # the process occurs as follows:
  # (_0) create a data.frame of relevant information and p.values
  # (_1) create col 'sig_P' with symbols for significance of P values
  # (_2) coerce col with symbols for significance of P values to character
  # (_3) round P values to the specified number of digits
  # (_4) reorder and retain relevant columns
  # (_5) rename columns 'P' and 'sig_P' using the input values for names_P
  if (inherits(data, "htest")) {
    kw_0 <- data.frame("Test" = "Kruskal-Wallis", "Comparison" = "global",
                       "P" = data$p.value, stringsAsFactors = F)
    kw_1 <- dplyr::mutate(kw_0, sig_P = symnum(P, cutpoints = cutpoints,
                                               symbols = symbols, corr = F))
    kw_2 <- dplyr::mutate(kw_1, sig_P = as.character(sig_P))
    kw_3 <- dplyr::mutate(kw_2, P = round(P, digits = digits))
    kw_4 <- dplyr::select(kw_3, Test, Comparison, P, sig_P)
    kw_5 <- dplyr::rename_at(kw_4, dplyr::vars(P, sig_P), ~names_P)
    return(kw_5)
  }
  # NOTE: the process is identical for both 'pairwise.htest' and 'PMCMR' ...
  # and occurs as follows:
  # (_0) create a data.frame of the p.values
  # (_1) create a data.frame of additional relevant information and p.values
  # (_2) reshape the data.frame into a workable format
  # (_3) remove any rows with NA vals (these represent same group comparisons)
  # (_4) rename column 'variable' to 'Pair2'
  # (_5) create col 'Comparison' by combining rows in cols 'Pair1' & 'Pair2'
  # (_6) create col 'sig_P' with symbols for significance of P values
  # (_7) coerce col with symbols for significance of P values to character
  # (_8) round P values to the specified number of digits
  # (_9) reorder and retain relevant columns
  # (_10) rename columns 'P' and 'sig_P' using the input values for names_P
  
  # input data are from pairwise.wilcox.test()
  if (inherits(data, "pairwise.htest")) {
    wi_0 <- data.frame(data$p.value, stringsAsFactors = F)
    wi_1 <- data.frame("Test" = "Wilcoxon", "Pair1" = row.names(wi_0), wi_0,
                       stringsAsFactors = F)
    wi_2 <- reshape2::melt(wi_1, id.vars = c("Test", "Pair1"), value.name = "P")
    wi_3 <- dplyr::filter(wi_2, !is.na(P))
    wi_4 <- dplyr::rename(wi_3, Pair2 = variable)
    wi_5 <- dplyr::mutate(wi_4,
                          Comparison = paste(Pair1, "vs", Pair2, sep = " "))
    wi_6 <- dplyr::mutate(wi_5, sig_P = symnum(P, cutpoints = cutpoints,
                                               symbols = symbols, corr = F))
    wi_7 <- dplyr::mutate(wi_6, sig_P = as.character(sig_P))
    wi_8 <- dplyr::mutate(wi_7, P = round(P, digits = digits))
    wi_9 <- dplyr::select(wi_8, Test, Comparison, P, sig_P)
    wi_10 <- dplyr::rename_at(wi_9, dplyr::vars(P, sig_P), ~names_P)
    return(wi_10)
  }
  # input data are from PMCMR::posthoc.kruskal.dunn.test()
  if (inherits(data, "PMCMR")) {
    du_0 <- data.frame(data$p.value, stringsAsFactors = F)
    du_1 <- data.frame("Test" = "Dunn's Test", "Pair1" = row.names(du_0), du_0,
                       stringsAsFactors = F)
    du_2 <- reshape2::melt(du_1, id.vars = c("Test", "Pair1"), value.name = "P")
    du_3 <- dplyr::filter(du_2, !is.na(P))
    du_4 <- dplyr::rename(du_3, Pair2 = variable)
    du_5 <- dplyr::mutate(du_4,
                          Comparison = paste(Pair1, "vs", Pair2, sep = " "))
    du_6 <- dplyr::mutate(du_5, sig_P = symnum(P, cutpoints = cutpoints,
                                               symbols = symbols, corr = F))
    du_7 <- dplyr::mutate(du_6, sig_P = as.character(sig_P))
    du_8 <- dplyr::mutate(du_7, P = round(P, digits = digits))
    du_9 <- dplyr::select(du_8, Test, Comparison, P, sig_P)
    du_10 <- dplyr::rename_at(du_9, dplyr::vars(P, sig_P), ~names_P)
    return(du_10)
  }
}
#
# example usage:
# new.df <- format_wi_kw(data = results, digits = 3,
#                        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
#                        symbols = c("****", "***", "**", "*", " ")
#                        names_P = c("BH_P", "sig_BH_P"))

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### E - STEP 1a - format data for statistics and plotting ----
### ************************************

# NOTE: if the environment is empty; there are some requirements:
# section_E requires access to objects in PREFACE; COSMOS; Section S

# provide provenance for information gathering at end of STEP 1i:
E.prov_secstep_ES1a <- "Section E - STEP 1a"
E.prov_heading_ES1a <- "format data for statistics and plotting"
E.prov_output_obj_ES1a <- "E.smp_adiv" # this object is output to the vault
E.prov_info <- paste("groups isolated from", E.prov_output_obj_ES1a, sep = " ")

# create new versions of the objects needed from Section S
# ... and create a vector naming those objects (used when saving the workspace)
E.sID_aom_EPD <- S.sID_aom_EPD
E.sID_aom_PD <- S.sID_aom_PD
E.sID_aom_E <- S.sID_aom_E
E.obj_from_S <- c("S.sID_aom_EPD",
                  "S.sID_aom_PD",
                  "S.sID_aom_E")

# format data to add sampledata to combined alpha diversity tables
# this process occurs as follows:
# (1) read in alpha diversity tables
# (2) rename column 'X' to SampleID
# (3:4) merge data.frames together and then merge with sampledata

E.raw_chao <- read.table(E.ifv_chao, header = T, sep = "\t", as.is = T,
                         stringsAsFactors = F)
E.raw_shan <- read.table(E.ifv_shan, header = T, sep = "\t", as.is = T,
                         stringsAsFactors = F)

E.chao <- dplyr::rename(E.raw_chao, SampleID = X)
E.shan <- dplyr::rename(E.raw_shan, SampleID = X)

E.adiv <- merge(x = E.chao, y = E.shan, sort = F, by = "SampleID")
E.smp_adiv <- merge(x = E.adiv, y = smp_dat, sort = F, by = "SampleID")

## use the SampleID vectors from Section S to subset the above data.frame
# for subset used for statistics, the process occurs as follows:
# (1:2) isolate samples and retain columns of interest
# (3) order by column DietType
# (4) convert column used as grouping variable into factors
# (5) convert to class data.frame
# (6) track provenance

# aom_EPD = AOM_DSS group cEcum; Proximal colon; Distal colon
E.adiv_aom_EPD_0 <- dplyr::filter(E.smp_adiv, SampleID %in% E.sID_aom_EPD)
E.adiv_aom_EPD_1 <- dplyr::select(E.adiv_aom_EPD_0, SampleID, chao1, shannon,
                                  Diet.Type, DietGroup, Type)
E.adiv_aom_EPD_2 <- E.adiv_aom_EPD_1[order(E.adiv_aom_EPD_1$Diet.Type), ]
E.adiv_aom_EPD_3 <- dplyr::mutate(
  E.adiv_aom_EPD_2, Diet.Type = factor(Diet.Type, levels = unique(Diet.Type)))
E.adiv_aom_EPD <- as.data.frame(E.adiv_aom_EPD_3)
E.info_adiv_aom_EPD <- data.frame(
  "object" = "E.adiv_aom_EPD",
  "value" = "AOM_DSS group cEcum; Proximal colon; Distal colon",
  stringsAsFactors = F)

# aom_cln = AOM_DSS group Proximal and Distal colon (grouped together)
E.adiv_aom_cln_0 <- dplyr::filter(E.smp_adiv, SampleID %in% E.sID_aom_PD)
E.adiv_aom_cln_1 <- dplyr::select(E.adiv_aom_cln_0, SampleID, chao1, shannon,
                                  DietGroup, DietGroupLab, Type, Diet.Type)
E.adiv_aom_cln_2 <- dplyr::mutate(
  E.adiv_aom_cln_1,
  DietGroupLab = factor(DietGroupLab, levels = unique(DietGroupLab)))
E.adiv_aom_cln <- E.adiv_aom_cln_2[order(E.adiv_aom_cln_2$DietGroup), ]
E.info_adiv_aom_cln <- data.frame(
  "object" = "E.adiv_aom_cln",
  "value" = "AOM_DSS group Proximal and Distal colon (grouped together)",
  stringsAsFactors = F)

# for subsets that will be plotted, the process occurs as follows:
# (1) define levels to order the specific groups plotted by ggplot2/ggpubr
# (2:3) isolate samples and retain the columns of interest
# (4) convert column used as grouping variable into factors with specific order
# (5) track provenance
# rinse and repeat all of the above steps for the other relevant subset(s)

# aom_E = AOM_DSS group cEcum
E.adiv_aom_diet_lvl_E <- c("C", "R", "F")
E.adiv_aom_E_0 <- dplyr::filter(E.smp_adiv, SampleID %in% E.sID_aom_E)
E.adiv_aom_E_1 <- dplyr::select(E.adiv_aom_E_0, SampleID, chao1, shannon,
                                DietGroup, DietGroupLab, Type)
E.adiv_aom_E <- dplyr::mutate(
  E.adiv_aom_E_1,
  DietGroupLab = factor(DietGroupLab, levels = E.adiv_aom_diet_lvl_E))
E.info_adiv_aom_E <- data.frame(
  "object" = "E.adiv_aom_E",
  "value" = "AOM_DSS group cEcum",
  stringsAsFactors = F)

# aom_PD = AOM_DSS group Proximal colon; Distal colon
E.adiv_aom_diet_lvl_PD <- c("C", "R", "F")
E.adiv_aom_dtyp_lvl_PD <- c("C.P", "C.D", "R.P", "R.D", "F.P", "F.D")
E.adiv_aom_PD_0 <- dplyr::filter(E.smp_adiv, SampleID %in% E.sID_aom_PD)
E.adiv_aom_PD_1 <- dplyr::select(E.adiv_aom_PD_0, SampleID, chao1, shannon,
                                 DietGroup, DietGroupLab, Type, Diet.Type)
E.adiv_aom_PD <- dplyr::mutate(
  E.adiv_aom_PD_1,
  DietGroupLab = factor(DietGroupLab, levels = E.adiv_aom_diet_lvl_PD), 
  Diet.Type = factor(Diet.Type, levels = E.adiv_aom_dtyp_lvl_PD))

E.info_adiv_aom_PD <- data.frame(
  "object" = "E.adiv_aom_PD",
  "value" = "AOM_DSS group Proximal colon; Distal colon",
  stringsAsFactors = F)

### ************************************
### E - STEP 1i - information gathering and provenance ----
### ************************************

# ** note for KDP: SS information gathering and provenance version 0.1 ** #

# (1) bind E.info data.frames together
# (2) add info/script/section/heading columns
# (3) reorder columns

E.info_0 <- rbind(E.info_adiv_aom_EPD,
                  E.info_adiv_aom_cln,
                  E.info_adiv_aom_PD,
                  E.info_adiv_aom_E)

E.info_1 <- E.info_0
E.info_1$info <- E.prov_info
E.info_1$script <- name_scrpt
E.info_1$section <- E.prov_secstep_ES1a
E.info_1$heading <- E.prov_heading_ES1a

E.info <- dplyr::select(E.info_1, info, value, object, script, section, heading)

### ************************************
### E - STEP 2a - statistics: perform global tests and format results ----
### ************************************

# perform 'global' test (i.e. All vs. All); this process occurs as follows:
# (1) perform test (NOTE: kruskal-wallis used as there are three groups)
# (2) pass test outputs to the format_wi_kw() function to format results
# (3:4) add new cols 'BH_P' and 'sig_BH_P' for merging w/ pwse results later on
# rinse and repeat (1:4) for the remaining relevant data.frame

# chao1
E.glbl_chao_aom_EPD_0 <- kruskal.test(chao1 ~ Diet.Type,
                                      data = E.adiv_aom_EPD)

E.glbl_chao_aom_EPD_1 <- format_wi_kw(data = E.glbl_chao_aom_EPD_0,
                                      digits = 5, cutpoints = cutpts,
                                      symbols = symbls,
                                      names_P = c("P", "sig_P"))

E.glbl_chao_aom_EPD <- E.glbl_chao_aom_EPD_1
E.glbl_chao_aom_EPD$BH_P <- "-"
E.glbl_chao_aom_EPD$sig_BH_P <- "-"

E.glbl_chao_aom_cln_0 <- kruskal.test(chao1 ~ DietGroupLab,
                                      data = E.adiv_aom_cln)

E.glbl_chao_aom_cln_1 <- format_wi_kw(data = E.glbl_chao_aom_cln_0,
                                      digits = 5, cutpoints = cutpts,
                                      symbols = symbls,
                                      names_P = c("P", "sig_P"))

E.glbl_chao_aom_cln <- E.glbl_chao_aom_cln_1
E.glbl_chao_aom_cln$BH_P <- "-"
E.glbl_chao_aom_cln$sig_BH_P <- "-"

# shannon
E.glbl_shan_aom_EPD_0 <- kruskal.test(shannon ~ Diet.Type,
                                      data = E.adiv_aom_EPD)

E.glbl_shan_aom_EPD_1 <- format_wi_kw(data = E.glbl_shan_aom_EPD_0,
                                      digits = 5, cutpoints = cutpts,
                                      symbols = symbls,
                                      names_P = c("P", "sig_P"))

E.glbl_shan_aom_EPD <- E.glbl_shan_aom_EPD_1
E.glbl_shan_aom_EPD$BH_P <- "-"
E.glbl_shan_aom_EPD$sig_BH_P <- "-"

E.glbl_shan_aom_cln_0 <- kruskal.test(shannon ~ DietGroupLab,
                                      data = E.adiv_aom_cln)

E.glbl_shan_aom_cln_1 <- format_wi_kw(data = E.glbl_shan_aom_cln_0,
                                      digits = 5, cutpoints = cutpts,
                                      symbols = symbls,
                                      names_P = c("P", "sig_P"))

E.glbl_shan_aom_cln <- E.glbl_shan_aom_cln_1
E.glbl_shan_aom_cln$BH_P <- "-"
E.glbl_shan_aom_cln$sig_BH_P <- "-"

### ************************************
### E - STEP 2b - statistics: perform pairwise tests and format results ----
### ************************************

# perform pairwise comparisons; this process occurs as follows:
# (1:2) perform pairwise tests without/with multiple comparisons adjustment
# (3) pass test outputs to the format_wi_kw() function to format results
# (4) merge appropriate _raw and _fdr data.frames by col 'Comparison'
# (5) combine appropriate _pwwi and _pwkw data.frames
# rinse and repeat (1:5) for the remaining relevant data.frame

# chao1
E.pwkw_raw_chao_aom_EPD_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_EPD$chao1,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "none")
E.pwwi_raw_chao_aom_EPD_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_EPD$chao1,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "none")
E.pwkw_fdr_chao_aom_EPD_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_EPD$chao1,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "fdr")
E.pwwi_fdr_chao_aom_EPD_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_EPD$chao1,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "fdr")

E.pwkw_raw_chao_aom_EPD_1 <- format_wi_kw(data = E.pwkw_raw_chao_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwkw_fdr_chao_aom_EPD_1 <- format_wi_kw(data = E.pwkw_fdr_chao_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))
E.pwwi_raw_chao_aom_EPD_1 <- format_wi_kw(data = E.pwwi_raw_chao_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwwi_fdr_chao_aom_EPD_1 <- format_wi_kw(data = E.pwwi_fdr_chao_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))

E.pwkw_chao_aom_EPD <- merge(x = E.pwkw_raw_chao_aom_EPD_1,
                             y = E.pwkw_fdr_chao_aom_EPD_1,
                             sort = F, by = c("Test", "Comparison"))
E.pwwi_chao_aom_EPD <- merge(x = E.pwwi_raw_chao_aom_EPD_1,
                             y = E.pwwi_fdr_chao_aom_EPD_1,
                             sort = F, by = c("Test", "Comparison"))

E.pwse_chao_aom_EPD <- rbind(E.pwkw_chao_aom_EPD, E.pwwi_chao_aom_EPD)

E.pwkw_raw_chao_aom_cln_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_cln$chao1,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "none")
E.pwwi_raw_chao_aom_cln_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_cln$chao1,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "none")
E.pwkw_fdr_chao_aom_cln_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_cln$chao1,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "fdr")
E.pwwi_fdr_chao_aom_cln_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_cln$chao1,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "fdr")

E.pwkw_raw_chao_aom_cln_1 <- format_wi_kw(data = E.pwkw_raw_chao_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwkw_fdr_chao_aom_cln_1 <- format_wi_kw(data = E.pwkw_fdr_chao_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))
E.pwwi_raw_chao_aom_cln_1 <- format_wi_kw(data = E.pwwi_raw_chao_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwwi_fdr_chao_aom_cln_1 <- format_wi_kw(data = E.pwwi_fdr_chao_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))

E.pwkw_chao_aom_cln <- merge(x = E.pwkw_raw_chao_aom_cln_1,
                             y = E.pwkw_fdr_chao_aom_cln_1,
                             sort = F, by = c("Test", "Comparison"))
E.pwwi_chao_aom_cln <- merge(x = E.pwwi_raw_chao_aom_cln_1,
                             y = E.pwwi_fdr_chao_aom_cln_1,
                             sort = F, by = c("Test", "Comparison"))

E.pwse_chao_aom_cln <- rbind(E.pwkw_chao_aom_cln, E.pwwi_chao_aom_cln)

# shannon
E.pwkw_raw_shan_aom_EPD_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_EPD$shannon,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "none")
E.pwwi_raw_shan_aom_EPD_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_EPD$shannon,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "none")
E.pwkw_fdr_shan_aom_EPD_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_EPD$shannon,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "fdr")
E.pwwi_fdr_shan_aom_EPD_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_EPD$shannon,
  g = E.adiv_aom_EPD$Diet.Type,
  p.adjust.method = "fdr")

E.pwkw_raw_shan_aom_EPD_1 <- format_wi_kw(data = E.pwkw_raw_shan_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwkw_fdr_shan_aom_EPD_1 <- format_wi_kw(data = E.pwkw_fdr_shan_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))
E.pwwi_raw_shan_aom_EPD_1 <- format_wi_kw(data = E.pwwi_raw_shan_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwwi_fdr_shan_aom_EPD_1 <- format_wi_kw(data = E.pwwi_fdr_shan_aom_EPD_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))

E.pwkw_shan_aom_EPD <- merge(x = E.pwkw_raw_shan_aom_EPD_1,
                             y = E.pwkw_fdr_shan_aom_EPD_1,
                             sort = F, by = c("Test", "Comparison"))
E.pwwi_shan_aom_EPD <- merge(x = E.pwwi_raw_shan_aom_EPD_1,
                             y = E.pwwi_fdr_shan_aom_EPD_1,
                             sort = F, by = c("Test", "Comparison"))

E.pwse_shan_aom_EPD <- rbind(E.pwkw_shan_aom_EPD, E.pwwi_shan_aom_EPD)


E.pwkw_raw_shan_aom_cln_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_cln$shannon,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "none")
E.pwwi_raw_shan_aom_cln_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_cln$shannon,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "none")
E.pwkw_fdr_shan_aom_cln_0 <- PMCMR::posthoc.kruskal.dunn.test(
  x = E.adiv_aom_cln$shannon,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "fdr")
E.pwwi_fdr_shan_aom_cln_0 <- pairwise.wilcox.test(
  x = E.adiv_aom_cln$shannon,
  g = E.adiv_aom_cln$DietGroupLab,
  p.adjust.method = "fdr")

E.pwkw_raw_shan_aom_cln_1 <- format_wi_kw(data = E.pwkw_raw_shan_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwkw_fdr_shan_aom_cln_1 <- format_wi_kw(data = E.pwkw_fdr_shan_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))
E.pwwi_raw_shan_aom_cln_1 <- format_wi_kw(data = E.pwwi_raw_shan_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("P", "sig_P"))
E.pwwi_fdr_shan_aom_cln_1 <- format_wi_kw(data = E.pwwi_fdr_shan_aom_cln_0,
                                          digits = 5, cutpoints = cutpts,
                                          symbols = symbls,
                                          names_P = c("BH_P", "sig_BH_P"))

E.pwkw_shan_aom_cln <- merge(x = E.pwkw_raw_shan_aom_cln_1,
                             y = E.pwkw_fdr_shan_aom_cln_1,
                             sort = F, by = c("Test", "Comparison"))
E.pwwi_shan_aom_cln <- merge(x = E.pwwi_raw_shan_aom_cln_1,
                             y = E.pwwi_fdr_shan_aom_cln_1,
                             sort = F, by = c("Test", "Comparison"))

E.pwse_shan_aom_cln <- rbind(E.pwkw_shan_aom_cln, E.pwwi_shan_aom_cln)

### ************************************
### E - STEP 2c - statistics: combine global and pairwise results ----
### ************************************

# provide provenance for information gathering at end of STEP 1i:
E.prov_secstep_ES2c <- "Section E - STEP 2c"
E.prov_heading_ES2c <- "statistics: combine global and pairwise results "
E.prov_outobj1_ES2c <- "E.rst_adiv_aom_EPD" # this object is output to the vault
E.prov_outobj2_ES2c <- "E.rst_flt_adiv_aom" # this object is output to the vault

# combine results data.frames as follows:
# (1) combine global and pwise results by metric
# (2;3) copy and add new column specifying metric
# (4) combine results together

E.rst_chao_aom_EPD_0 <- rbind(E.glbl_chao_aom_EPD, E.pwse_chao_aom_EPD)
E.rst_shan_aom_EPD_0 <- rbind(E.glbl_shan_aom_EPD, E.pwse_shan_aom_EPD)
E.rst_chao_aom_cln_0 <- rbind(E.glbl_chao_aom_cln, E.pwse_chao_aom_cln)
E.rst_shan_aom_cln_0 <- rbind(E.glbl_shan_aom_cln, E.pwse_shan_aom_cln)
E.rst_chao_aom_EPD <- E.rst_chao_aom_EPD_0
E.rst_shan_aom_EPD <- E.rst_shan_aom_EPD_0
E.rst_chao_aom_cln <- E.rst_chao_aom_cln_0
E.rst_shan_aom_cln <- E.rst_shan_aom_cln_0
E.rst_chao_aom_EPD$metric <- "chao1"
E.rst_shan_aom_EPD$metric <- "shannon"
E.rst_chao_aom_cln$metric <- "chao1"
E.rst_shan_aom_cln$metric <- "shannon"
E.rst_adiv_aom_EPD_0 <- rbind(E.rst_chao_aom_EPD, E.rst_shan_aom_EPD)
E.rst_adiv_aom_cln_0 <- rbind(E.rst_chao_aom_cln, E.rst_shan_aom_cln)

# for _EPD, format the data as follows:
# (1) define a vector to isolate between group comparisons (btw)
# (2) define a vector to isolate within group comparisons (wtn)
# (3) isolate the comparisons; NOTE: (xtr) = extras
# (4:5) copy and add new column to specify category (i.e. between/within/extra)
# (6) reorder columns
# (7) combine results back together
# (8) order by relevant columns

E.btw_comp <- c("R.E vs C.E", "F.E vs C.E", "R.E vs F.E",
                "R.P vs C.P", "F.P vs C.P", "R.P vs F.P",
                "R.D vs C.D", "F.D vs C.D", "R.D vs F.D")

E.wtn_comp <- c("C.P vs C.E", "C.E vs C.D", "C.P vs C.D",
                "R.P vs R.E", "R.E vs R.D", "R.P vs R.D",
                "F.P vs F.E", "F.E vs F.D", "F.P vs F.D")

E.btw_rst_adiv_aom_EPD_0 <- dplyr::filter(E.rst_adiv_aom_EPD_0,
                                          Comparison %in% E.btw_comp)
E.wtn_rst_adiv_aom_EPD_0 <- dplyr::filter(E.rst_adiv_aom_EPD_0,
                                          Comparison %in% E.wtn_comp)
E.xtr_rst_adiv_aom_EPD_0 <- dplyr::filter(E.rst_adiv_aom_EPD_0,
                                          !Comparison %in% E.btw_comp,
                                          !Comparison %in% E.wtn_comp)

E.btw_rst_adiv_aom_EPD_1 <- E.btw_rst_adiv_aom_EPD_0
E.wtn_rst_adiv_aom_EPD_1 <- E.wtn_rst_adiv_aom_EPD_0
E.xtr_rst_adiv_aom_EPD_1 <- E.xtr_rst_adiv_aom_EPD_0
E.btw_rst_adiv_aom_EPD_1$category <- "between"
E.wtn_rst_adiv_aom_EPD_1$category <- "within"
E.xtr_rst_adiv_aom_EPD_1$category <- "extra"

E.btw_rst_adiv_aom_EPD <- dplyr::select(E.btw_rst_adiv_aom_EPD_1,
                                        metric, category, Comparison,
                                        dplyr::everything())
E.wtn_rst_adiv_aom_EPD <- dplyr::select(E.wtn_rst_adiv_aom_EPD_1,
                                        metric, category, Comparison,
                                        dplyr::everything())
E.xtr_rst_adiv_aom_EPD <- dplyr::select(E.xtr_rst_adiv_aom_EPD_1,
                                        metric, category, Comparison,
                                        dplyr::everything())

E.rst_adiv_aom_EPD_1 <- rbind(E.btw_rst_adiv_aom_EPD, E.wtn_rst_adiv_aom_EPD,
                              E.xtr_rst_adiv_aom_EPD)
E.rst_adiv_aom_EPD <- E.rst_adiv_aom_EPD_1[
  order(E.rst_adiv_aom_EPD_1$metric,
        E.rst_adiv_aom_EPD_1$Test,
        E.rst_adiv_aom_EPD_1$category,
        E.rst_adiv_aom_EPD_1$Comparison), ]

## filter results
# for _EPD, this process occurs as follows:
# (1) format to retain only the btw and wtn comparisons for only Wilcoxon
# (2:3) copy & and coerce column Comparison into a character
# (4:6) copy & add full descriptions for info Diet.Type info in col Comparison

E.rfl_adiv_aom_EPD_0 <- dplyr::filter(E.rst_adiv_aom_EPD,
                                      Test == "Wilcoxon",
                                      category == "between" |
                                        category == "within")
E.rfl_adiv_aom_EPD_1 <- E.rfl_adiv_aom_EPD_0
E.rfl_adiv_aom_EPD_1$Comparison <- as.character(E.rfl_adiv_aom_EPD_1$Comparison)
E.rfl_adiv_aom_EPD <- E.rfl_adiv_aom_EPD_1
E.rfl_adiv_aom_EPD$CompKey <- E.rfl_adiv_aom_EPD$Comparison
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = ".E",
                                   replacement = ".cecum",
                                   x = E.rfl_adiv_aom_EPD$CompKey)
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = ".P",
                                   replacement = ".proximal colon",
                                   x = E.rfl_adiv_aom_EPD$CompKey)
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = ".D",
                                   replacement = ".distal colon",
                                   x = E.rfl_adiv_aom_EPD$CompKey)
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = "C.",
                                   replacement = "Control ",
                                   x = E.rfl_adiv_aom_EPD$CompKey)
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = "R.",
                                   replacement = "Rice bran ",
                                   x = E.rfl_adiv_aom_EPD$CompKey)
E.rfl_adiv_aom_EPD$CompKey <- gsub(pattern = "F.",
                                   replacement = "Fermented rice bran ",
                                   x = E.rfl_adiv_aom_EPD$CompKey)

# for _cln, this process occurs as follows:
# (1:2) copy and add column category
# (3) format to retain only Wilcoxon
# (4) copy & add full descriptions for DietGroup info in column Comparison

E.rst_adiv_aom_cln_1 <- E.rst_adiv_aom_cln_0
E.rst_adiv_aom_cln_1$category <- "between"
E.rst_adiv_aom_cln_2 <- dplyr::filter(E.rst_adiv_aom_cln_1,
                                      Test == "Wilcoxon")
E.rfl_adiv_aom_cln <- E.rst_adiv_aom_cln_2
E.rfl_adiv_aom_cln$CompKey <- gsub(
  pattern = "C vs F",
  replacement = "Control versus Fermented rice bran (colon combined)",
  x = E.rfl_adiv_aom_cln$Comparison)
E.rfl_adiv_aom_cln$CompKey <- gsub(
  pattern = "R vs C",
  replacement = "Rice bran versus Control (colon combined)",
  x = E.rfl_adiv_aom_cln$CompKey)
E.rfl_adiv_aom_cln$CompKey <- gsub(
  pattern = "R vs F",
  replacement = "Rice bran versus Fermented rice bran (colon combined)",
  x = E.rfl_adiv_aom_cln$CompKey)

# combine the _rfl data.frames and reorder columns
E.rst_flt_adiv_aom_0 <- rbind(E.rfl_adiv_aom_cln, E.rfl_adiv_aom_EPD)
E.rst_flt_adiv_aom <- dplyr::select(E.rst_flt_adiv_aom_0,
                                    metric, category, Comparison, Test,
                                    P, sig_P, BH_P, sig_BH_P, CompKey)

### ************************************
### E - STEP  3 - format for the addition of BH-P values/symbols to plots ----
### ************************************

# NOTE: the stat_compare_means() function will be used to plot brackets...
# ... for each comparison while labels will originate from the BH-P values ...
# ... generated in STEP 2 above (i.e. use info in E.rst_flt_adiv_aom)
# NOTE: proceeding this way also allows us to plot the BH/fdr adjusted sigs ...
# ... rather than the undajusted values used by stat_compare_means()

# define comparisons for statistical tests that are used in plots
E.comps_grp <- list(c("C", "R"), c("R", "F"), c("C", "F"))

# define vectors to isolate:
# (1) between group comparisons for cecum (cec)
# (2) between group comparisons for colon (cln)

E.cec_comp <- c("R.E vs C.E", "R.E vs F.E", "F.E vs C.E")
E.cln_comp <- c("R vs C", "R vs F", "C vs F")

# isolate the above comparins to create data.frames to add BH-P values
E.cec_bhp_chao_0 <- dplyr::filter(E.rst_flt_adiv_aom, metric == "chao1",
                                  Comparison %in% E.cec_comp)
E.cec_bhp_shan_0 <- dplyr::filter(E.rst_flt_adiv_aom, metric == "shannon",
                                  Comparison %in% E.cec_comp)
E.cln_bhp_chao_0 <- dplyr::filter(E.rst_flt_adiv_aom, metric == "chao1",
                                  Comparison %in% E.cln_comp)
E.cln_bhp_shan_0 <- dplyr::filter(E.rst_flt_adiv_aom, metric == "shannon",
                                  Comparison %in% E.cln_comp)

# define scalers to scale y value for labels (i.e. places labels above brackets)
E.scl_chao <- 18
E.scl_shan <- 0.36

# (1:2) define x and y locations for each comparison
# (3:) copy and add the x and y locations for each comparison

E.x_CvR <- 1.5
E.x_RvF <- 2.5
E.x_CvF <- 2.0
E.y_chao_CvR <- 360
E.y_chao_RvF <- 400
E.y_chao_CvF <- 440
E.y_shan_CvR <- 7.2
E.y_shan_RvF <- 8.0
E.y_shan_CvF <- 8.8

E.cec_bhp_chao_1 <- E.cec_bhp_chao_0
E.cec_bhp_chao_1$x[E.cec_bhp_chao_1$Comparison == "R.E vs C.E"] <- E.x_CvR
E.cec_bhp_chao_1$x[E.cec_bhp_chao_1$Comparison == "R.E vs F.E"] <- E.x_RvF
E.cec_bhp_chao_1$x[E.cec_bhp_chao_1$Comparison == "F.E vs C.E"] <- E.x_CvF
E.cec_bhp_chao_1$y[E.cec_bhp_chao_1$Comparison == "R.E vs C.E"] <- E.y_chao_CvR
E.cec_bhp_chao_1$y[E.cec_bhp_chao_1$Comparison == "R.E vs F.E"] <- E.y_chao_RvF
E.cec_bhp_chao_1$y[E.cec_bhp_chao_1$Comparison == "F.E vs C.E"] <- E.y_chao_CvF

E.cec_bhp_shan_1 <- E.cec_bhp_shan_0
E.cec_bhp_shan_1$x[E.cec_bhp_shan_1$Comparison == "R.E vs C.E"] <- E.x_CvR
E.cec_bhp_shan_1$x[E.cec_bhp_shan_1$Comparison == "R.E vs F.E"] <- E.x_RvF
E.cec_bhp_shan_1$x[E.cec_bhp_shan_1$Comparison == "F.E vs C.E"] <- E.x_CvF
E.cec_bhp_shan_1$y[E.cec_bhp_shan_1$Comparison == "R.E vs C.E"] <- E.y_shan_CvR
E.cec_bhp_shan_1$y[E.cec_bhp_shan_1$Comparison == "R.E vs F.E"] <- E.y_shan_RvF
E.cec_bhp_shan_1$y[E.cec_bhp_shan_1$Comparison == "F.E vs C.E"] <- E.y_shan_CvF

E.cln_bhp_chao_1 <- E.cln_bhp_chao_0
E.cln_bhp_chao_1$x[E.cln_bhp_chao_1$Comparison == "R vs C"] <- E.x_CvR
E.cln_bhp_chao_1$x[E.cln_bhp_chao_1$Comparison == "R vs F"] <- E.x_RvF
E.cln_bhp_chao_1$x[E.cln_bhp_chao_1$Comparison == "C vs F"] <- E.x_CvF
E.cln_bhp_chao_1$y[E.cln_bhp_chao_1$Comparison == "R vs C"] <- E.y_chao_CvR
E.cln_bhp_chao_1$y[E.cln_bhp_chao_1$Comparison == "R vs F"] <- E.y_chao_RvF
E.cln_bhp_chao_1$y[E.cln_bhp_chao_1$Comparison == "C vs F"] <- E.y_chao_CvF

E.cln_bhp_shan_1 <- E.cln_bhp_shan_0
E.cln_bhp_shan_1$x[E.cln_bhp_shan_1$Comparison == "R vs C"] <- E.x_CvR
E.cln_bhp_shan_1$x[E.cln_bhp_shan_1$Comparison == "R vs F"] <- E.x_RvF
E.cln_bhp_shan_1$x[E.cln_bhp_shan_1$Comparison == "C vs F"] <- E.x_CvF
E.cln_bhp_shan_1$y[E.cln_bhp_shan_1$Comparison == "R vs C"] <- E.y_shan_CvR
E.cln_bhp_shan_1$y[E.cln_bhp_shan_1$Comparison == "R vs F"] <- E.y_shan_RvF
E.cln_bhp_shan_1$y[E.cln_bhp_shan_1$Comparison == "C vs F"] <- E.y_shan_CvF

# (1) define a vector to label the BH_P value
# (2:4) copy & coerce BH_P to numeric and then run a logical test to create lab

E.bhp_lab <- "BH-P = "

E.cec_bhp_chao_2 <- E.cec_bhp_chao_1
E.cec_bhp_shan_2 <- E.cec_bhp_shan_1
E.cln_bhp_chao_2 <- E.cln_bhp_chao_1
E.cln_bhp_shan_2 <- E.cln_bhp_shan_1
E.cec_bhp_chao_2$BH_P <- as.numeric(E.cec_bhp_chao_2$BH_P)
E.cec_bhp_shan_2$BH_P <- as.numeric(E.cec_bhp_shan_2$BH_P)
E.cln_bhp_chao_2$BH_P <- as.numeric(E.cln_bhp_chao_2$BH_P)
E.cln_bhp_shan_2$BH_P <- as.numeric(E.cln_bhp_shan_2$BH_P)

E.cec_bhp_chao_2$lab_BH_P <- ifelse(E.cec_bhp_chao_2$BH_P < 0.05, 
                                    yes = paste(sep = "", E.bhp_lab, 
                                                E.cec_bhp_chao_2$BH_P), 
                                    no = "ns")
E.cec_bhp_shan_2$lab_BH_P <- ifelse(E.cec_bhp_shan_2$BH_P < 0.05, 
                                    yes = paste(sep = "", E.bhp_lab, 
                                                E.cec_bhp_shan_2$BH_P), 
                                    no = "ns")
E.cln_bhp_chao_2$lab_BH_P <- ifelse(E.cln_bhp_chao_2$BH_P < 0.05, 
                                    yes = paste(sep = "", E.bhp_lab, 
                                                E.cln_bhp_chao_2$BH_P), 
                                    no = "ns")
E.cln_bhp_shan_2$lab_BH_P <- ifelse(E.cln_bhp_shan_2$BH_P < 0.05, 
                                    yes = paste(sep = "", E.bhp_lab, 
                                                E.cln_bhp_shan_2$BH_P), 
                                    no = "ns")

# (1) create a new column to label with the symbol
E.cec_bhp_chao_3 <- dplyr::mutate(
  E.cec_bhp_chao_2, 
  sym_BH_P = symnum(BH_P, 
                    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), 
                    corr = F))

E.cec_bhp_shan_3 <- dplyr::mutate(
  E.cec_bhp_shan_2, 
  sym_BH_P = symnum(BH_P, 
                    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), 
                    corr = F))

E.cln_bhp_chao_3 <- dplyr::mutate(
  E.cln_bhp_chao_2, 
  sym_BH_P = symnum(BH_P, 
                    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), 
                    corr = F))

E.cln_bhp_shan_3 <- dplyr::mutate(
  E.cln_bhp_shan_2, 
  sym_BH_P = symnum(BH_P, 
                    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), 
                    corr = F))

# copy to rename and coerce col 'sym_BH_P' to character
E.cec_bhp_chao_4 <- E.cec_bhp_chao_3
E.cec_bhp_shan_4 <- E.cec_bhp_shan_3
E.cln_bhp_chao_4 <- E.cln_bhp_chao_3
E.cln_bhp_shan_4 <- E.cln_bhp_shan_3
E.cec_bhp_chao_4$sym_BH_P <- as.character(E.cec_bhp_chao_4$sym_BH_P)
E.cec_bhp_shan_4$sym_BH_P <- as.character(E.cec_bhp_shan_4$sym_BH_P)
E.cln_bhp_chao_4$sym_BH_P <- as.character(E.cln_bhp_chao_4$sym_BH_P)
E.cln_bhp_shan_4$sym_BH_P <- as.character(E.cln_bhp_shan_4$sym_BH_P)

# copy to rename final
E.cec_bhp_chao <- E.cec_bhp_chao_4
E.cec_bhp_shan <- E.cec_bhp_shan_4
E.cln_bhp_chao <- E.cln_bhp_chao_4
E.cln_bhp_shan <- E.cln_bhp_shan_4

### ************************************
### E - STEP  4 - define plot parameters/customize plot aesthetics ----
### ************************************

# NOTE: for DietGroupLab aesthetics, order = C, R, F
# NOTE: for Diet.Type aesthetics, order = C.P, C.D, R.P, R.D, F.P, F.D
# these orders are defined at the end of "E - STEP 1a"

# define x axis breaks and labels
E.xbrk <- c("C", "R", "F")
E.xlab <- c(expression(~Con^spf), expression(~RB^spf), expression(~FRB^spf))

# define y axis limits, breaks, and labels
E.ylim_chao <- c(0, 460)
E.ybrk_chao <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)
E.ylab_chao <- c("0", "", " 100", "", " 200", "", " 300", "", " 400")

E.ylim_shan <- c(0, 9.4)
E.ybrk_shan <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
E.ylab_shan <- c("0", "", "2", "", "4", "", "6", "", "8")

## shapes, colors, labels
# (1) vector for cecum DietGroupLab shapes
# (2) vector for cecum DietGroupLab colors
# (3) vector for cecum DietGroupLab fill
# (4) vector for colon Diet.Type shapes
# (5) vector for colon Diet.Type colors
# (6) vector for colon Diet.Type fill
# (7) vector for DietGroupLab box outline colors

E.shp_diet_cec <- c(shp_C, shp_R, shp_F)
E.hex_diet_cec <- c(hex_C[1], hex_R[1], hex_F[1])
E.fil_diet_cec <- c(hex_C[2], hex_R[2], hex_F[2])
E.shp_diet_cln <- c(shp_C, shp_C, shp_R, shp_R, shp_F, shp_F)
E.hex_diet_cln <- c(hex_C[1], hex_C[2], hex_R[1], hex_R[2], hex_F[1], hex_F[2])
E.fil_diet_cln <- c(hex_C[2], hex_C[4], hex_R[2], hex_R[4], hex_F[2], hex_F[4])
E.bxp_diet <- c(hex_C[3], hex_R[3], hex_F[3])

## define jitter w/ set random seed to ensure reproducibility of final plots
# (1) jitter for cecum
# (2) jitter and dodge for colon

E.jitr_cec <- position_jitter(width = 0.13, seed = 123456789)
E.jitr_cln <- position_jitterdodge(jitter.width = 0.42, dodge.width = 0.51,
                                   seed = 123456789)

## sizing parameters
E.sze_sts_lab <- 3.00 # size for BH-P text labels
E.sze_tip_len <- 0.00 # size for bracket tip length
E.sze_bkt_lne <- 0.50 # size for bracket lines
E.sze_smp_pts <- 1.51 # size for sample points
E.sze_bxp_wid <- 0.55 # size for width of boxes
E.sze_bxp_lwd <- 0.36 # size for box lines
E.sze_ptl <- 16 # size for plot title (i.e. Type)
E.sze_stl <- 14 # size for plot title (i.e. metric)
E.sze_atx <- 12 # size for axis text (x and y)

# plot labeling
E.ttl_cec <- text_grob("Cecum", face = "bold", color = greydient[1], 
                       size = E.sze_ptl)
E.ttl_cln <- text_grob("Colon (proximal and distal)", face = "bold",
                       color = greydient[1], size = E.sze_ptl)
E.ytl_chao <- "Chao1 Index"
E.ytl_shan <- "Shannon Index"

# custom ggplot theme parameters
E.JittrBox <- theme(
  axis.ticks = element_line(color = greydient[1]),
  axis.text = element_text(size = E.sze_atx),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = E.sze_atx, hjust = 0.5),
  plot.title = element_text(size = E.sze_stl, hjust = 0.5),
  legend.position = "none",
  panel.ontop = T,
  panel.background = element_rect(fill = NA),
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]))

### ************************************
### E - STEP 5a - plotting: create stripcharts with boxplots ----
### ************************************

# provide provenance for information gathering at end of STEP 1i:
E.prov_secstep_ES5a <- "Section E - STEP 5a"
E.prov_heading_ES5a <- "plotting: stripcharts with boxplots"
# NOTE: each of the plots below are output to the vault

# begin plotting
E.gpr_chao_cec <- ggstripchart(data = E.adiv_aom_E, ylab = E.ytl_chao, 
                               x = "DietGroupLab", y = "chao1",
                               shape = "DietGroupLab",
                               fill = "DietGroupLab",
                               color = "DietGroupLab",
                               size = E.sze_smp_pts,
                               font.family = fnt, 
                               position = E.jitr_cec, order = E.xbrk, 
                               add = "boxplot", 
                               add.params = list(width = E.sze_bxp_wid,
                                                 size = E.sze_bxp_lwd,
                                                 alpha = 0,
                                                 color = E.bxp_diet,
                                                 fill = greydient[8])) +
  stat_compare_means(inherit.aes = T, comparison = E.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = E.sze_bkt_lne, tip.length = E.sze_tip_len,
                     label.y = c(E.y_chao_CvR, E.y_chao_RvF, E.y_chao_CvF)) +
  # number label
  # geom_text(data = E.cec_bhp_chao, family = fnt, color = greydient[2],
  #           size = E.sze_sts_lab, inherit.aes = F,
  #           aes(x = x, y = y + E.scl_chao, label = lab_BH_P)) +
  # symbol label
  geom_text(data = E.cec_bhp_chao, family = fnt, color = greydient[2],
            size = E.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + E.scl_chao, label = sym_BH_P)) +
  scale_shape_manual(values = E.shp_diet_cec) +
  scale_color_manual(values = E.hex_diet_cec) +
  scale_fill_manual(values = E.fil_diet_cec) +
  scale_x_discrete(breaks = E.xbrk, labels = E.xlab) +
  scale_y_continuous(limits = E.ylim_chao, breaks = E.ybrk_chao, 
                     labels = E.ylab_chao) +
  E.JittrBox 

E.gpr_shan_cec <- ggstripchart(data = E.adiv_aom_E, ylab = E.ytl_shan,
                               x = "DietGroupLab", y = "shannon",
                               shape = "DietGroupLab",
                               fill = "DietGroupLab",
                               color = "DietGroupLab",
                               size = E.sze_smp_pts,
                               font.family = fnt, 
                               position = E.jitr_cec, order = E.xbrk, 
                               add = "boxplot", 
                               add.params = list(width = E.sze_bxp_wid,
                                                 size = E.sze_bxp_lwd,
                                                 alpha = 0,
                                                 color = E.bxp_diet,
                                                 fill = greydient[8])) +
  stat_compare_means(inherit.aes = T, comparison = E.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = E.sze_bkt_lne, tip.length = E.sze_tip_len,
                     label.y = c(E.y_shan_CvR, E.y_shan_RvF, E.y_shan_CvF)) +
  geom_text(data = E.cec_bhp_shan, family = fnt, color = greydient[2],
            size = E.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + E.scl_shan, label = sym_BH_P)) +
  scale_shape_manual(values = E.shp_diet_cec) +
  scale_color_manual(values = E.hex_diet_cec) +
  scale_fill_manual(values = E.fil_diet_cec) +
  scale_x_discrete(breaks = E.xbrk, labels = E.xlab) +
  scale_y_continuous(limits = E.ylim_shan, breaks = E.ybrk_shan, 
                     labels = E.ylab_shan) +
  E.JittrBox

E.gpr_chao_cln <- ggstripchart(data = E.adiv_aom_PD, ylab = E.ytl_chao,
                               x = "DietGroupLab", y = "chao1",
                               shape = "Diet.Type",
                               fill = "Diet.Type",
                               color = "Diet.Type",
                               size = E.sze_smp_pts,
                               font.family = fnt, 
                               position = E.jitr_cln, order = E.xbrk,
                               add = "boxplot", 
                               add.params = list(width = E.sze_bxp_wid,
                                                 size = E.sze_bxp_lwd,
                                                 alpha = 0,
                                                 color = E.bxp_diet,
                                                 fill = greydient[8])) +
  stat_compare_means(inherit.aes = T, comparison = E.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = E.sze_bkt_lne, tip.length = E.sze_tip_len,
                     label.y = c(E.y_chao_CvR, E.y_chao_RvF, E.y_chao_CvF)) +
  geom_text(data = E.cln_bhp_chao, family = fnt, color = greydient[2],
            size = E.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + E.scl_chao, label = sym_BH_P)) +
  scale_shape_manual(values = E.shp_diet_cln) +
  scale_color_manual(values = E.hex_diet_cln) +
  scale_fill_manual(values = E.fil_diet_cln) +
  scale_x_discrete(breaks = E.xbrk, labels = E.xlab) +
  scale_y_continuous(limits = E.ylim_chao, breaks = E.ybrk_chao, 
                     labels = E.ylab_chao) +
  E.JittrBox

E.gpr_shan_cln <- ggstripchart(data = E.adiv_aom_PD, ylab = E.ytl_shan, 
                               x = "DietGroupLab", y = "shannon",
                               shape = "Diet.Type",
                               fill = "Diet.Type",
                               color = "Diet.Type",
                               size = E.sze_smp_pts,
                               font.family = fnt, 
                               position = E.jitr_cln, order = E.xbrk,
                               add = "boxplot", 
                               add.params = list(width = E.sze_bxp_wid,
                                                 size = E.sze_bxp_lwd,
                                                 alpha = 0,
                                                 color = E.bxp_diet,
                                                 fill = greydient[8])) +
  stat_compare_means(inherit.aes = T, comparison = E.comps_grp, geom = "text",
                     method = "wilcox.test", size = 0, 
                     bracket.size = E.sze_bkt_lne, tip.length = E.sze_tip_len,
                     label.y = c(E.y_shan_CvR, E.y_shan_RvF, E.y_shan_CvF)) +
  geom_text(data = E.cln_bhp_shan, family = fnt, color = greydient[2],
            size = E.sze_sts_lab, inherit.aes = F,
            aes(x = x, y = y + E.scl_shan, label = sym_BH_P)) +
  scale_shape_manual(values = E.shp_diet_cln) +
  scale_color_manual(values = E.hex_diet_cln) +
  scale_fill_manual(values = E.fil_diet_cln) +
  scale_x_discrete(breaks = E.xbrk, labels = E.xlab) +
  scale_y_continuous(limits = E.ylim_shan, breaks = E.ybrk_shan, 
                     labels = E.ylab_shan) +
  E.JittrBox

### ************************************
### E - STEP 5b - plotting: arrange plots and format ----
### ************************************

# provide provenance for information gathering at end of STEP 1i:
E.prov_secstep_ES5b <- "Section E - STEP 5b"
E.prov_heading_ES5b <- "plotting: arrange plots and format"
# NOTE: each of the plots below are output to the vault

# the arrangement process occurs as follows:
# (1) arrange by Type
# (2) define Type labels
# (4) annotate arrangements with Type labels and panel labels
# rinse and repeat

E.gga_cec_0 <- ggarrange(E.gpr_chao_cec, E.gpr_shan_cec, labels = NULL,
                         ncol = 2, nrow = 1, align = "hv")
E.gga_cec <- annotate_figure(p = E.gga_cec_0, top = E.ttl_cec,
                             fig.lab = "A", fig.lab.face = "bold",
                             fig.lab.size = 30)

E.gga_cln_0 <- ggarrange(E.gpr_chao_cln, E.gpr_shan_cln, labels = NULL,
                         ncol = 2, nrow = 1, align = "hv")
E.gga_cln <- annotate_figure(p = E.gga_cln_0, top = E.ttl_cln,
                             fig.lab = "B", fig.lab.face = "bold",
                             fig.lab.size = 30)

# # arrange formatted plots together
# E.gga_adiv <- ggarrange(E.gga_cec, E.gga_cln, labels = c("",""),
#                         ncol = 2, nrow = 1, align = "hv")
# ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
#        filename = "plot.pdf", plot = E.gga_adiv)

### ************************************
### E - WRITE OUTPUTS ----
### ************************************

# provenance for outputs to the vault
E.prov_output1 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_smp_adiv,
                             "object" = E.prov_output_obj_ES1a,
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES1a,
                             "heading" = E.prov_heading_ES1a,
                             stringsAsFactors = F)
E.prov_output2 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_rst_adiv_EPD,
                             "object" = E.prov_outobj1_ES2c,
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES2c,
                             "heading" = E.prov_heading_ES2c,
                             stringsAsFactors = F)
E.prov_output3 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_rst_flt_adiv,
                             "object" = E.prov_outobj2_ES2c,
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES2c,
                             "heading" = E.prov_heading_ES2c,
                             stringsAsFactors = F)
E.prov_output4 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gpr_chao_cec,
                             "object" = "E.gpr_chao_C_SNTE",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5a,
                             "heading" = E.prov_heading_ES5a,
                             stringsAsFactors = F)
E.prov_output5 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gpr_shan_cec,
                             "object" = "E.gpr_shan_C_SNTE",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5a,
                             "heading" = E.prov_heading_ES5a,
                             stringsAsFactors = F)
E.prov_output6 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gpr_chao_cln,
                             "object" = "E.gpr_chao_R_SNTE",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5a,
                             "heading" = E.prov_heading_ES5a,
                             stringsAsFactors = F)
E.prov_output7 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gpr_shan_cln,
                             "object" = "E.gpr_shan_R_SNTE",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5a,
                             "heading" = E.prov_heading_ES5a,
                             stringsAsFactors = F)
E.prov_output8 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gga_cec,
                             "object" = "E.gga_cec",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5b,
                             "heading" = E.prov_heading_ES5b,
                             stringsAsFactors = F)
E.prov_output9 <- data.frame("info" = "provenance for output",
                             "path" = E.ofv_gga_cln,
                             "object" = "E.gga_cln",
                             "script" = name_scrpt,
                             "section" = E.prov_secstep_ES5b,
                             "heading" = E.prov_heading_ES5b,
                             stringsAsFactors = F)
E.prov <- rbind(E.prov_output1, E.prov_output2, E.prov_output3, E.prov_output4, 
                E.prov_output5, E.prov_output6, E.prov_output7, E.prov_output8,
                E.prov_output9)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = E.smp_adiv, file = E.ofv_smp_adiv)
write.table(sep = "\t", row.names = F, x = E.rst_adiv_aom_EPD, 
            file = E.ofv_rst_adiv_EPD)
write.table(sep = "\t", row.names = F, x = E.rst_flt_adiv_aom,
            file = E.ofv_rst_flt_adiv)

write.table(sep = "\t", row.names = F, x = E.info, file = E.ofv_info)
write.table(sep = "\t", row.names = F, x = E.prov, file = E.ofv_prov)

ggsave(device = "pdf", dpi = 600, units = "mm", width = 82.5, height = 90,
       filename = E.ofv_gpr_chao_cec, plot = E.gpr_chao_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 82.5, height = 90,
       filename = E.ofv_gpr_shan_cec, plot = E.gpr_shan_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 82.5, height = 90,
       filename = E.ofv_gpr_chao_cln, plot = E.gpr_chao_cln)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 82.5, height = 90,
       filename = E.ofv_gpr_shan_cln, plot = E.gpr_shan_cln)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = E.ofv_gga_cec, plot = E.gga_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = E.ofv_gga_cln, plot = E.gga_cln)

E.obj <- ls(pattern = "E.")
E.lst <- c(E.obj[grep(pattern = "E.", x = E.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS, E.obj_from_S)
save(list = E.lst, file = E.ofv_wksp)

### ************************************
### B - temporary preface ----
### ************************************

# setwd("~/Desktop/") # *gasp*
# load(file = "SarumansScroll/vault/WS_Script02_PREFACE.RData") # load PREFACE
# load(file = S.ofv_wksp) # load Section S
# load(file = ofv_COSMOS_wksp) # load COSMOS
# require(ggplot2)
# require(ggpubr)
T.EBTKS_pro_tree <- ape::read.tree(T.ofv_EBTKS_pro_tree) # read in Sec. T's tree

### ************************************
### B - FUNCTIONS ----
### ************************************

# create a vector naming all Section B functions (used when saving workspace)
B.functions <- c("iso_var_exp", "iso_cord", "flt_smp_dst", "pwise_adon")

# ** note for KDP: version 0.2 ** #
# function takes an input of class "pcoa" or class "prcomp" and isolates...
# the amount of variation explained by the PC specified
# NOTE: isolates 1 value at a time as specified by pc.num
# NOTE: returns % values with 3 digits
# iso_var_exp = isolate_variation_explained
iso_var_exp <- function(data = c(pcoa, prcomp), pc.num, 
                        output = c("value", "label")) {
  # internal check to ensure correct input class
  if (!inherits(data, c("pcoa", "prcomp"))) {
    stop("data must be class 'pcoa' OR class 'prcomp'")
  }
  if (!inherits(pc.num, "numeric")) {
    stop("pc.num must be numeric")
  }
  if (!output == "value" && !output == "label") {
    stop("output must be one of: 'value' OR 'label'")
  }
  if (inherits(data, "pcoa")) {
    # isolate pcoa$values (create a data.frame)
    pcoa_vals <- data.frame(data$values)
    # isolate the Relative_eig values
    rel_eigs <- pcoa_vals$Relative_eig
    # convert rel eigs to percentage and round to 3 digits
    rel_eigs_rnd <- (round(rel_eigs, digits = 3)) * 100
    # isolate the specified PC
    pcoa_val <- rel_eigs_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pcoa_val)
    }
    if (output == "label") {
      pcoa_lab <- paste("PC ", pc.num, " (", pcoa_val, "%)", sep = "")
      return(pcoa_lab)
    }
  }
  if (inherits(data, "prcomp")) {
    # isolate standard deviations for all pcs
    pcs_sdv <- data$sdev
    # square each sdev, sum the squares, and divide each sdev by the sum
    pcs_sdv_sqr <- pcs_sdv^2
    pcs_sdv_sqr_sum <- sum(pcs_sdv_sqr)
    pca_vals_raw <- pcs_sdv_sqr/pcs_sdv_sqr_sum
    # convert to percentage and round to 3 digits
    pca_vals_rnd <- (round(pca_vals_raw, digits = 3)) * 100
    # isolate the specified PC
    pca_val <- pca_vals_rnd[pc.num]
    # statements to determine what to return
    if (output == "value") {
      return(pca_val)
    }
    if (output == "label") {
      pca_lab <- paste("PC ", pc.num, " (", pca_val, "%)", sep = "")
      return(pca_lab)
    }
  }
}
# 
# # example usage:
# val <- iso_var_exp(data = pcoa, pc.num = 1, output = "value") # returns number
# lab <- iso_var_exp(data = pcoa, pc.num = 1, output = "label") # returns label

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.2 ** #
# function takes an input of class "pcoa" isolates...
# the x and y coordinates for each sample
# iso_cord = isolate_coordinates
iso_cord <- function(data = pcoa) {
  # internal check to ensure correct input class
  if (!inherits(data, "pcoa")) {
    stop("data must be class 'pcoa'")
  }
  # isolate pcoa$vectors (create a data.frame)
  pcoa_vec <- data.frame(data$vectors)
  # isolate cols 1 and 2 (Axis.1, Axis.2; i.e. x and y coordinates)
  pcoa_crd <- pcoa_vec[, 1:2]
  # rename cols
  pcoa_crd_rnm <- pcoa_crd
  names(pcoa_crd_rnm)[1] <- "xvar" # version 0.2 changed "X.Coord" to "xvar"
  names(pcoa_crd_rnm)[2] <- "yvar" # version 0.2 changed "Y.Coord" to "yvar"
  return(pcoa_crd_rnm)
}
#
# example usage:
# pcoa_cord_df <- iso_cord(data = pcoa)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.1 ** #
# function takes inputs of sampledata and a distance matrix and...
# ... returns filtered versions of both output as a list
# the resultant outputs contain identical samples which is useful when...
# ... running PERMANOVA/adonis or using the pwise_adon() function
# NOTE: this function is a dependency for pwise_adon()
flt_smp_dst <- function(smp_dat = data.frame, dst_obj = dist, flt_col = "", 
                        flt_row = "", match_col = "") {
  # internal checks to ensure correct input classes
  if (!inherits(smp_dat, "data.frame")) {
    stop("input smp_dat must be class 'data.frame'")
  }
  if (!inherits(dst_obj, "dist")) {
    stop("input dst_obj must be class 'dist'")
  }
  if (!inherits(flt_col, "character")) {
    stop("flt_col must be character")
  }
  if (!inherits(flt_row, "character")) {
    stop("flt_row must be character")
  }
  if (!inherits(match_col, "character")) {
    stop("match_col must be character")
  }
  # create vector of row #s in sampledata flt_col that match the flt_row
  row_num <- which(smp_dat[, flt_col] %in% flt_row)
  if (isTRUE(length(row_num) == 0)) {
    stop("`flt_row` not found in `flt_col`")
  }
  # create new sampledata df using the above vector
  new_smp_dat_0 <- smp_dat[row_num, ]
  # coerce distance matrix from class: 'dist' to class: 'matrix'
  dmx <- as.matrix(dst_obj)
  # create a vector of row #s in the dmx that match column in the sampledata
  dmx_vec_nums <- which(row.names(dmx) %in% new_smp_dat_0[, match_col])
  # filter the distance matrix to retain only the numbers in the above vector
  # NOTE: we can take advantage of matrixland as col #s and row #s are identical
  new_dmx <- dmx[dmx_vec_nums, dmx_vec_nums]
  # match the order rows in new sampledata and new distance matrix
  new_smp_dat_1 <- new_smp_dat_0[match(x = row.names(new_dmx), 
                                       table = new_smp_dat_0[, match_col]), ]
  # convert new_dmx back into class 'dist'
  new_dst <- as.dist(new_dmx)
  # create a list containing new_smp_dat and new.dist
  new_lst <- list(new_smp = new_smp_dat_1, new_dst = new_dst)
  #return(new_lst)
  return(new_lst)
}
#
# example usage:
# filtered_list <- flt_smp_dst(smp_dat = smp, dst_obj = dst, flt_col = "group",
#                              flt_row = "housecat", match_col = "SampleID")
# filtered_smp_dat <- filtered_list$new_smp
# filtered_dst_obj <- filtered_list$new_dst

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# ** note for KDP: version 0.3 ** #
# function takes inputs of sampledata and a distance matrix and...
# ... performs pairwise PERMANOVA comparisons for all possible pairs...
# ... present in the specified input 'test_col' and returns a data.frame with...
# ... the comparison, R2 value, unadjusted p-values & p values adjusted...
# ... using the specified adjustment method
# NOTE: this function is dependent on the flt_smp_dst() function
pwise_adon <- function(smp_dat = data.frame, dst_obj = dist, test_col = "", 
                       match_col = "", perms, digits,
                       adjust = c("holm", "hochberg", "hommel", "bonferroni", 
                                  "BH", "BY", "fdr", "none")) {
  # internal checks to ensure correct input classes
  if (!inherits(smp_dat, "data.frame")) {
    stop("input smp_dat must be class 'data.frame'")
  }
  if (!inherits(dst_obj, "dist")) {
    stop("input dst_obj must be class 'dist'")
  }
  if (!inherits(test_col, "character")) {
    stop("test_col must be character")
  }
  if (!inherits(match_col, "character")) {
    stop("match_col must be character")
  }
  if (!adjust == "holm" && !adjust == "hochberg" && !adjust == "hommel" &&
      !adjust == "bonferroni" && !adjust == "BH" && !adjust == "BY" &&
      !adjust == "fdr" && !adjust == "none") {
    stop("adjust must be one of the methods supported by p.adjust();
    'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'")
  }
  # create a vector of row values in sampledata flt_col
  row_val <- smp_dat[, test_col]
  # create unique factors for row values in sampledata flt_col
  factors <- as.factor(unique(row_val))
  # obtain all pairwise comparisons from the above row values
  comps <- combn(x = unique(row_val), m = 2)
  # create empty vectors to store values obtained in the for loop below
  pair_labs <- c()
  R2 <- c()
  p_raw <- c()
  # loop through the data to make adonis comparisons for each pair in comps
  for(i in 1:ncol(comps)){
    # create a vector for an individual pair
    pair <- as.vector(factors[factors %in% c(comps[1, i], comps[2, i])])
    # filter the sampledata df and the distance matrix for the paired comparison
    filt <- flt_smp_dst(flt_row = pair, smp_dat = smp_dat, dst_obj = dst_obj, 
                        flt_col = test_col, match_col = match_col)
    # create filtered data
    flt_smp <- filt$new_smp
    # create new distance matrix reflecting the filtering step above
    flt_dst <- filt$new_dst
    # perform adonis on the filtered data
    adon <- vegan::adonis(formula = flt_dst ~ flt_smp[, test_col], 
                          data = flt_smp, permutations = perms)
    # define labels for the individual pair
    pair_labs <- c(pair_labs, paste(pair[1], 'vs', pair[2]))
    # store information from the adonis data.frame
    R2 <- c(R2, adon$aov.tab[1, 5])
    p_raw <- c(p_raw, adon$aov.tab[1, 6])
  }
  # adjust raw p values for multiple comparisons
  # NOTE: number of comparisons (n) is determined by length(p_raw)
  p_adj <- p.adjust(p = p_raw, method = adjust, n = length(p_raw))
  # create new data.frame containing the relvant information
  pwise_df_0 <- data.frame("Comparison" = pair_labs, R2, p_raw, p_adj)
  # append input adjustment method to column name for p_adj
  pwise_df_1 <- pwise_df_0
  p_adj_num <- which(names(pwise_df_1) == "p_adj")
  names(pwise_df_1)[p_adj_num] <- paste("p_", adjust, sep = "")
  # round all numeric values to specified number of decimal places
  pwise_df_rnd <- dplyr::mutate_if(pwise_df_1, is.numeric, 
                                   list(~ round(., digits)))
  return(pwise_df_rnd)
}
#
# example usage:
# results_pwise_adon <- pwise_adon(smp_dat = smp, dst_obj = dst,
#                                  test_col = "group", match_col = "SampleID",
#                                  perms = 999, adjust = "BH", digits = 3)

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

### ************************************
### B - STEP  1 - format data for UniFracs/CZM/clr/PCA ----
### ************************************

# NOTE: if the environment is empty; there are some requirements:
# section_B requires access to objects in PREFACE; COSMOS; Section S; Section T

# create new versions of the objects needed from sections S and T...
# ... and create a vector naming those objects (used when saving the workspace)
B.EBTKS_pro_tree <- T.EBTKS_pro_tree
B.abs_aom_EPD_0 <- S.abs_aom_EPD
B.abs_aom_PD_0 <- S.abs_aom_PD
B.abs_aom_E_0 <- S.abs_aom_E
B.obj_from_ST <- c("T.EBTKS_pro_tree",
                   "S.abs_aom_EPD", "S.abs_aom_PD", "S.abs_aom_E")

# NOTE: aom_EPD is for statistics; aom_E is for plotting; aom_PD is for both

# format data to obtain RepSeq as columns and SampleIDs as rows
# this process occurs as follows:
# (1) copy to avoid overwriting data
# (2) convert column RepSeq into row.names
# (3) remove unneeded columns
# (4) transpose to convert features to cols and samples to rows
# rinse and repeat (1:4) for other objects

B.abs_aom_EPD_1 <- B.abs_aom_EPD_0
row.names(B.abs_aom_EPD_1) <- B.abs_aom_EPD_1$RepSeq
B.abs_aom_EPD_2 <- dplyr::select(B.abs_aom_EPD_1, -dplyr::one_of(com_col))
B.abs_aom_EPD <- as.data.frame(t(B.abs_aom_EPD_2))

B.abs_aom_PD_1 <- B.abs_aom_PD_0
row.names(B.abs_aom_PD_1) <- B.abs_aom_PD_1$RepSeq
B.abs_aom_PD_2 <- dplyr::select(B.abs_aom_PD_1, -dplyr::one_of(com_col))
B.abs_aom_PD <- as.data.frame(t(B.abs_aom_PD_2))

B.abs_aom_E_1 <- B.abs_aom_E_0
row.names(B.abs_aom_E_1) <- B.abs_aom_E_1$RepSeq
B.abs_aom_E_2 <- dplyr::select(B.abs_aom_E_1, -dplyr::one_of(com_col))
B.abs_aom_E <- as.data.frame(t(B.abs_aom_E_2))

### ************************************
### B - STEP  2 - compute UniFracs/PCoA ----
### ************************************

# for _EPD and _cln, the process occurs as follows:
# (1) compute UniFracs (NOTE: input data = features as cols and samples as rows)
# (2) extract relevant metrics from the unifrac lists

B.uni_aom_EPD <- GUniFrac::GUniFrac(otu.tab = B.abs_aom_EPD, 
                                    tree = B.EBTKS_pro_tree,
                                    alpha = c(0, 0.5, 1))$unifracs
B.guni_aom_EPD <- as.dist(B.uni_aom_EPD[, , "d_0.5"])
B.uuni_aom_EPD <- as.dist(B.uni_aom_EPD[, , "d_UW"])
B.wuni_aom_EPD <- as.dist(B.uni_aom_EPD[, , "d_1"])

B.uni_aom_cln <- GUniFrac::GUniFrac(otu.tab = B.abs_aom_PD, 
                                    tree = B.EBTKS_pro_tree,
                                    alpha = c(0, 0.5, 1))$unifracs
B.guni_aom_cln <- as.dist(B.uni_aom_cln[, , "d_0.5"])
B.uuni_aom_cln <- as.dist(B.uni_aom_cln[, , "d_UW"])
B.wuni_aom_cln <- as.dist(B.uni_aom_cln[, , "d_1"])

# for subsets that will be plotted, the process occurs as follows:
# (1) compute UniFracs (NOTE: input data = features as cols and samples as rows)
# (2) extract Generalized UniFrac from the unifrac lists
# (3) compute PCoA
# (4) create PCoA plot axis labels (i.e. isolate variation explained by PC1;PC2)
# (5) isolate x & y coordinates
# (6:8) convert row.names (i.e. SampleIDs) into a col & merge with sampledata
# (9) retain columns of interest
# (10) define levels to order the specific groups plotted by ggplot2/ggpubr
# (11) convert column used as grouping variable into factors with specific order
# rinse and repeat all of the above steps for the other relevant subset(s)

B.uni_aom_E <- GUniFrac::GUniFrac(otu.tab = B.abs_aom_E, 
                                  tree = B.EBTKS_pro_tree,
                                  alpha = c(0, 0.5, 1))$unifracs
B.guni_aom_E <- as.dist(B.uni_aom_E[, , "d_0.5"])
B.pcoa_guni_aom_E <- ape::pcoa(B.guni_aom_E)
B.pc1_guni_aom_E <- iso_var_exp(data = B.pcoa_guni_aom_E, pc.num = 1,
                                output = "label")
B.pc2_guni_aom_E <- iso_var_exp(data = B.pcoa_guni_aom_E, pc.num = 2,
                                output = "label")
B.cord_guni_aom_E_0 <- iso_cord(data = B.pcoa_guni_aom_E)
B.cord_guni_aom_E_1 <- B.cord_guni_aom_E_0
B.cord_guni_aom_E_1$SampleID <- row.names(B.cord_guni_aom_E_1)
B.smp_cord_guni_aom_E_0 <- merge(x = B.cord_guni_aom_E_1, y = smp_dat,
                                 all = F, sort = F, by = "SampleID")
B.smp_cord_guni_aom_E_1 <- dplyr::select(B.smp_cord_guni_aom_E_0, 
                                         SampleID, xvar, yvar,
                                         DietGroupLab, Type)
B.guni_aom_diet_lvl_E <- c("C", "R", "F")
B.smp_cord_guni_aom_E <- dplyr::mutate(
  B.smp_cord_guni_aom_E_1,
  DietGroupLab = factor(DietGroupLab, levels = B.guni_aom_diet_lvl_E))

B.uni_aom_PD <- GUniFrac::GUniFrac(otu.tab = B.abs_aom_PD, 
                                   tree = B.EBTKS_pro_tree,
                                   alpha = c(0, 0.5, 1))$unifracs
B.guni_aom_PD <- as.dist(B.uni_aom_PD[, , "d_0.5"])
B.pcoa_guni_aom_PD <- ape::pcoa(B.guni_aom_PD)
B.pc1_guni_aom_PD <- iso_var_exp(data = B.pcoa_guni_aom_PD, pc.num = 1,
                                 output = "label")
B.pc2_guni_aom_PD <- iso_var_exp(data = B.pcoa_guni_aom_PD, pc.num = 2,
                                 output = "label")
B.cord_guni_aom_PD_0 <- iso_cord(data = B.pcoa_guni_aom_PD)
B.cord_guni_aom_PD_1 <- B.cord_guni_aom_PD_0
B.cord_guni_aom_PD_1$SampleID <- row.names(B.cord_guni_aom_PD_1)
B.smp_cord_guni_aom_PD_0 <- merge(x = B.cord_guni_aom_PD_1, y = smp_dat,
                                  all = F, sort = F, by = "SampleID")
B.smp_cord_guni_aom_PD_1 <- dplyr::select(B.smp_cord_guni_aom_PD_0, 
                                          SampleID, xvar, yvar,
                                          DietGroupLab, Type, Diet.Type)
B.guni_aom_diet_lvl_PD <- c("C", "R", "F")
B.guni_aom_dtyp_lvl_PD <- c("C.P", "C.D", "R.P", "R.D", "F.P", "F.D")
B.smp_cord_guni_aom_PD <- dplyr::mutate(
  B.smp_cord_guni_aom_PD_1,
  DietGroupLab = factor(DietGroupLab, levels = B.guni_aom_diet_lvl_PD), 
  Diet.Type = factor(Diet.Type, levels = B.guni_aom_dtyp_lvl_PD))

### ************************************
### B - STEP  3 - compute CZM/clr/PCA ----
### ************************************

# for CZM and clr, the process occurs as follows:
# (1) replace zero counts
# (2) transpose to convert features to rows and samples to cols
# (3) clr transform
# rinse and repeat (1:3) for other objects
# NOTE: for cmultRepl(), input data = features as cols and samples as rows
# NOTE: for clr transform, input data = features as rows and samples as cols

B.czm_aom_EPD_0 <- zCompositions::cmultRepl(B.abs_aom_EPD, method = "CZM", 
                                            output = "prop")
B.czm_aom_EPD_1 <- as.data.frame(t(B.czm_aom_EPD_0))
B.clr_aom_EPD_0 <- as.data.frame(
  apply(B.czm_aom_EPD_1, 2, function(x) {log2(x) - mean(log2(x))}))

B.czm_aom_cln_0 <- zCompositions::cmultRepl(B.abs_aom_PD, method = "CZM", 
                                            output = "prop")
B.czm_aom_cln_1 <- as.data.frame(t(B.czm_aom_cln_0))
B.clr_aom_cln_0 <- as.data.frame(
  apply(B.czm_aom_cln_1, 2, function(x) {log2(x) - mean(log2(x))}))

B.czm_aom_PD_0 <- zCompositions::cmultRepl(B.abs_aom_PD, method = "CZM", 
                                           output = "prop")
B.czm_aom_PD_1 <- as.data.frame(t(B.czm_aom_PD_0))
B.clr_aom_PD_0 <- as.data.frame(
  apply(B.czm_aom_PD_1, 2, function(x) {log2(x) - mean(log2(x))}))

B.czm_aom_E_0 <- zCompositions::cmultRepl(B.abs_aom_E, method = "CZM", 
                                          output = "prop")
B.czm_aom_E_1 <- as.data.frame(t(B.czm_aom_E_0))
B.clr_aom_E_0 <- as.data.frame(
  apply(B.czm_aom_E_1, 2, function(x) {log2(x) - mean(log2(x))}))

# for PCA, the process occurs as follows:
# (1) transpose to convert features to cols and samples to rows
# (2) compute PCA
# (3) create PCA plot axis labels (i.e. isolate variation explained by PC1;PC2)
# (4:6) create ggbiplots, extract x & y coordinates, merge with sampledata
# (7) retain columns of interest (and rename column 'labels' to 'SampleID')
# (8) order by column DietGroupLab and Type
# (9) define levels to order the specific groups plotted by ggplot2/ggpubr
# (10) convert column used as grouping variable into factors with specific order
# rinse and repeat all of the above steps for the other relevant subset
# NOTE: for prcomp(), input data = features as cols and samples as rows

B.clr_aom_E_1 <- as.data.frame(t(B.clr_aom_E_0))
B.pca_atch_aom_E <- prcomp(B.clr_aom_E_1)
B.pc1_atch_aom_E <- iso_var_exp(B.pca_atch_aom_E, pc.num = 1, 
                                output = "label")
B.pc2_atch_aom_E <- iso_var_exp(B.pca_atch_aom_E, pc.num = 2, 
                                output = "label")
B.gbp_pca_atch_aom_E <- ggbiplot::ggbiplot(
  B.pca_atch_aom_E, labels = row.names(B.clr_aom_E_1), scale = 0, 
  var.axes = F) # print(B.gbp_pca_atch_aom_E)
B.cord_atch_aom_E <- B.gbp_pca_atch_aom_E[["data"]]
B.smp_cord_atch_aom_E_0 <- merge(x = B.cord_atch_aom_E, y = smp_dat, all = F,
                                 sort = F, by.x = "labels", by.y = "SampleID")
B.smp_cord_atch_aom_E_1 <- dplyr::select(B.smp_cord_atch_aom_E_0,
                                         SampleID = labels, xvar, yvar,
                                         DietGroupLab, Type)
B.atch_aom_diet_lvl_E <- c("C", "R", "F")
B.smp_cord_atch_aom_E <- dplyr::mutate(
  B.smp_cord_atch_aom_E_1,
  DietGroupLab = factor(DietGroupLab, levels = B.atch_aom_diet_lvl_E))

B.clr_aom_PD_1 <- as.data.frame(t(B.clr_aom_PD_0))
B.pca_atch_aom_PD <- prcomp(B.clr_aom_PD_1)
B.pc1_atch_aom_PD <- iso_var_exp(B.pca_atch_aom_PD, pc.num = 1, 
                                 output = "label")
B.pc2_atch_aom_PD <- iso_var_exp(B.pca_atch_aom_PD, pc.num = 2, 
                                 output = "label")
B.gbp_pca_atch_aom_PD <- ggbiplot::ggbiplot(
  B.pca_atch_aom_PD, labels = row.names(B.clr_aom_PD_1), scale = 0, 
  var.axes = F) # print(B.gbp_pca_atch_aom_PD)
B.cord_atch_aom_PD <- B.gbp_pca_atch_aom_PD[["data"]]
B.smp_cord_atch_aom_PD_0 <- merge(x = B.cord_atch_aom_PD, y = smp_dat, all = F,
                                  sort = F, by.x = "labels", by.y = "SampleID")
B.smp_cord_atch_aom_PD_1 <- dplyr::select(B.smp_cord_atch_aom_PD_0, 
                                          SampleID = labels, xvar, yvar,
                                          DietGroupLab, Type, Diet.Type)
B.atch_aom_diet_lvl_PD <- c("C", "R", "F")
B.atch_aom_dtyp_lvl_PD <- c("C.P", "C.D", "R.P", "R.D", "F.P", "F.D")
B.smp_cord_atch_aom_PD <- dplyr::mutate(
  B.smp_cord_atch_aom_PD_1,
  DietGroupLab = factor(DietGroupLab, levels = B.atch_aom_diet_lvl_PD), 
  Diet.Type = factor(Diet.Type, levels = B.atch_aom_dtyp_lvl_PD))

### ************************************
### B - STEP 4a - statistics: perform pairwise PERMANOVA ----
### ************************************

# for CZM/clr dfs, the process occurs as follows:
# (1) transpose to convert features to cols and samples to rows
# (2) create distance matrices using the euclidean metric
# NOTE: for dist(), input data = features as cols and samples as rows

B.clr_aom_EPD_1 <- as.data.frame(t(B.clr_aom_EPD_0))
B.atch_aom_EPD <- dist(B.clr_aom_EPD_1, method = "euclidean")

B.clr_aom_cln_1 <- as.data.frame(t(B.clr_aom_cln_0))
B.atch_aom_cln <- dist(B.clr_aom_cln_1, method = "euclidean")

# filter sampledata & distance matrices to have identical samples
# NOTE: the goal here is just to create a filtered smp_dat object for each comp
B.list_guni_aom_EPD <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.guni_aom_EPD, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_uuni_aom_EPD <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uuni_aom_EPD, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_wuni_aom_EPD <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wuni_aom_EPD, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_atch_aom_EPD <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.atch_aom_EPD, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")

B.smp_guni_aom_EPD <- B.list_guni_aom_EPD$new_smp
B.dst_guni_aom_EPD <- B.list_guni_aom_EPD$new_dst
B.smp_uuni_aom_EPD <- B.list_uuni_aom_EPD$new_smp
B.dst_uuni_aom_EPD <- B.list_uuni_aom_EPD$new_dst
B.smp_wuni_aom_EPD <- B.list_wuni_aom_EPD$new_smp
B.dst_wuni_aom_EPD <- B.list_wuni_aom_EPD$new_dst
B.smp_atch_aom_EPD <- B.list_atch_aom_EPD$new_smp
B.dst_atch_aom_EPD <- B.list_atch_aom_EPD$new_dst

B.list_guni_aom_cln <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.guni_aom_cln, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_uuni_aom_cln <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.uuni_aom_cln, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_wuni_aom_cln <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.wuni_aom_cln, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")
B.list_atch_aom_cln <- flt_smp_dst(smp_dat = smp_dat, dst_obj = B.atch_aom_cln, 
                                   flt_col = "Cohort", flt_row = "AOM_DSS", 
                                   match_col = "SampleID")

B.smp_guni_aom_cln <- B.list_guni_aom_cln$new_smp
B.dst_guni_aom_cln <- B.list_guni_aom_cln$new_dst
B.smp_uuni_aom_cln <- B.list_uuni_aom_cln$new_smp
B.dst_uuni_aom_cln <- B.list_uuni_aom_cln$new_dst
B.smp_wuni_aom_cln <- B.list_wuni_aom_cln$new_smp
B.dst_wuni_aom_cln <- B.list_wuni_aom_cln$new_dst
B.smp_atch_aom_cln <- B.list_atch_aom_cln$new_smp
B.dst_atch_aom_cln <- B.list_atch_aom_cln$new_dst

# perform pairwise PERMANOVA using adonis
B.prm_guni_aom_EPD_0 <- pwise_adon(smp_dat = B.smp_guni_aom_EPD, 
                                   dst_obj = B.dst_guni_aom_EPD,
                                   test_col = "Diet.Type",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_uuni_aom_EPD_0 <- pwise_adon(smp_dat = B.smp_uuni_aom_EPD, 
                                   dst_obj = B.dst_uuni_aom_EPD,
                                   test_col = "Diet.Type",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_wuni_aom_EPD_0 <- pwise_adon(smp_dat = B.smp_wuni_aom_EPD, 
                                   dst_obj = B.dst_wuni_aom_EPD,
                                   test_col = "Diet.Type",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_atch_aom_EPD_0 <- pwise_adon(smp_dat = B.smp_atch_aom_EPD, 
                                   dst_obj = B.dst_atch_aom_EPD,
                                   test_col = "Diet.Type",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)

B.prm_guni_aom_cln_0 <- pwise_adon(smp_dat = B.smp_guni_aom_cln, 
                                   dst_obj = B.dst_guni_aom_cln,
                                   test_col = "DietGroupLab",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_uuni_aom_cln_0 <- pwise_adon(smp_dat = B.smp_uuni_aom_cln, 
                                   dst_obj = B.dst_uuni_aom_cln,
                                   test_col = "DietGroupLab",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_wuni_aom_cln_0 <- pwise_adon(smp_dat = B.smp_wuni_aom_cln, 
                                   dst_obj = B.dst_wuni_aom_cln,
                                   test_col = "DietGroupLab",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)
B.prm_atch_aom_cln_0 <- pwise_adon(smp_dat = B.smp_atch_aom_cln, 
                                   dst_obj = B.dst_atch_aom_cln,
                                   test_col = "DietGroupLab",
                                   match_col = "SampleID", perms = 9999, 
                                   adjust = "BH", digits = 4)

### ************************************
### B - STEP 4b - statistics: format PERMANOVA results ----
### ************************************

# provide provenance for information gathering at end of STEP 4i:
B.prov_secstep_BS4b <- "Section B - STEP 4b"
B.prov_heading_BS4b <- "statistics: format PERMANOVA results"
B.prov_outobj1_BS4b <- "B.rst_bdiv_aom_EPD" # this object is output to the vault
B.prov_outobj2_BS4b <- "B.rst_flt_bdiv_aom" # this object is output to the vault

# for _EPD and_cln, format the data as follows:
# (1:2) copy and add new column specifying metric

B.prm_guni_aom_EPD <- B.prm_guni_aom_EPD_0
B.prm_uuni_aom_EPD <- B.prm_uuni_aom_EPD_0
B.prm_wuni_aom_EPD <- B.prm_wuni_aom_EPD_0
B.prm_atch_aom_EPD <- B.prm_atch_aom_EPD_0
B.prm_guni_aom_cln <- B.prm_guni_aom_cln_0
B.prm_uuni_aom_cln <- B.prm_uuni_aom_cln_0
B.prm_wuni_aom_cln <- B.prm_wuni_aom_cln_0
B.prm_atch_aom_cln <- B.prm_atch_aom_cln_0

B.prm_guni_aom_EPD$metric <- "generalized UniFrac"
B.prm_uuni_aom_EPD$metric <- "unweighted UniFrac"
B.prm_wuni_aom_EPD$metric <- "weighted UniFrac"
B.prm_atch_aom_EPD$metric <- "aitchison/euclidean"
B.prm_guni_aom_cln$metric <- "generalized UniFrac"
B.prm_uuni_aom_cln$metric <- "unweighted UniFrac"
B.prm_wuni_aom_cln$metric <- "weighted UniFrac"
B.prm_atch_aom_cln$metric <- "aitchison/euclidean"

# for _EPD, format the data as follows:
# (1) define a vector to isolate between group comparisons (btw)
# (2) define a vector to isolate within group comparisons (wtn)
# (3) isolate the comparisons; NOTE: (xtr) = extras
# (4:6) define vector to be used as row values in new column category
# (7:8) copy and add new column to specify category (i.e. between/within/extra)
# (9) reorder columns
# (10:13) reorder columns, order data, and then bind btw with wtn and xtr
# (14) combine results back together

B.btw_comp <- c("C.E vs R.E", "C.E vs F.E", "R.E vs F.E",
                "C.P vs R.P", "C.P vs F.P", "R.P vs F.P",
                "C.D vs R.D", "C.D vs F.D", "R.D vs F.D")

B.wtn_comp <- c("C.P vs C.E", "C.D vs C.E", "C.P vs C.D",
                "R.P vs R.E", "R.D vs R.E", "R.P vs R.D",
                "F.P vs F.E", "F.D vs F.E", "F.P vs F.D")

B.btw_rst_guni_aom_EPD_0 <- dplyr::filter(B.prm_guni_aom_EPD,
                                          Comparison %in% B.btw_comp)
B.wtn_rst_guni_aom_EPD_0 <- dplyr::filter(B.prm_guni_aom_EPD,
                                          Comparison %in% B.wtn_comp)
B.xtr_rst_guni_aom_EPD_0 <- dplyr::filter(B.prm_guni_aom_EPD,
                                          !Comparison %in% B.btw_comp, 
                                          !Comparison %in% B.wtn_comp)
B.btw_rst_uuni_aom_EPD_0 <- dplyr::filter(B.prm_uuni_aom_EPD,
                                          Comparison %in% B.btw_comp)
B.wtn_rst_uuni_aom_EPD_0 <- dplyr::filter(B.prm_uuni_aom_EPD,
                                          Comparison %in% B.wtn_comp)
B.xtr_rst_uuni_aom_EPD_0 <- dplyr::filter(B.prm_uuni_aom_EPD,
                                          !Comparison %in% B.btw_comp, 
                                          !Comparison %in% B.wtn_comp)
B.btw_rst_wuni_aom_EPD_0 <- dplyr::filter(B.prm_wuni_aom_EPD,
                                          Comparison %in% B.btw_comp)
B.wtn_rst_wuni_aom_EPD_0 <- dplyr::filter(B.prm_wuni_aom_EPD,
                                          Comparison %in% B.wtn_comp)
B.xtr_rst_wuni_aom_EPD_0 <- dplyr::filter(B.prm_wuni_aom_EPD,
                                          !Comparison %in% B.btw_comp, 
                                          !Comparison %in% B.wtn_comp)
B.btw_rst_atch_aom_EPD_0 <- dplyr::filter(B.prm_atch_aom_EPD,
                                          Comparison %in% B.btw_comp)
B.wtn_rst_atch_aom_EPD_0 <- dplyr::filter(B.prm_atch_aom_EPD,
                                          Comparison %in% B.wtn_comp)
B.xtr_rst_atch_aom_EPD_0 <- dplyr::filter(B.prm_atch_aom_EPD,
                                          !Comparison %in% B.btw_comp, 
                                          !Comparison %in% B.wtn_comp)

B.btw_row <- "between"
B.wtn_row <- "within"
B.xtr_row <- "extra"

B.btw_rst_guni_aom_EPD_1 <- B.btw_rst_guni_aom_EPD_0
B.wtn_rst_guni_aom_EPD_1 <- B.wtn_rst_guni_aom_EPD_0
B.xtr_rst_guni_aom_EPD_1 <- B.xtr_rst_guni_aom_EPD_0
B.btw_rst_uuni_aom_EPD_1 <- B.btw_rst_uuni_aom_EPD_0
B.wtn_rst_uuni_aom_EPD_1 <- B.wtn_rst_uuni_aom_EPD_0
B.xtr_rst_uuni_aom_EPD_1 <- B.xtr_rst_uuni_aom_EPD_0
B.btw_rst_wuni_aom_EPD_1 <- B.btw_rst_wuni_aom_EPD_0
B.wtn_rst_wuni_aom_EPD_1 <- B.wtn_rst_wuni_aom_EPD_0
B.xtr_rst_wuni_aom_EPD_1 <- B.xtr_rst_wuni_aom_EPD_0
B.btw_rst_atch_aom_EPD_1 <- B.btw_rst_atch_aom_EPD_0
B.wtn_rst_atch_aom_EPD_1 <- B.wtn_rst_atch_aom_EPD_0
B.xtr_rst_atch_aom_EPD_1 <- B.xtr_rst_atch_aom_EPD_0
B.btw_rst_guni_aom_EPD_1$category <- B.btw_row
B.wtn_rst_guni_aom_EPD_1$category <- B.wtn_row
B.xtr_rst_guni_aom_EPD_1$category <- B.xtr_row
B.btw_rst_uuni_aom_EPD_1$category <- B.btw_row
B.wtn_rst_uuni_aom_EPD_1$category <- B.wtn_row
B.xtr_rst_uuni_aom_EPD_1$category <- B.xtr_row
B.btw_rst_wuni_aom_EPD_1$category <- B.btw_row
B.wtn_rst_wuni_aom_EPD_1$category <- B.wtn_row
B.xtr_rst_wuni_aom_EPD_1$category <- B.xtr_row
B.btw_rst_atch_aom_EPD_1$category <- B.btw_row
B.wtn_rst_atch_aom_EPD_1$category <- B.wtn_row
B.xtr_rst_atch_aom_EPD_1$category <- B.xtr_row

B.btw_rst_guni_aom_EPD_2 <- dplyr::select(B.btw_rst_guni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.wtn_rst_guni_aom_EPD_2 <- dplyr::select(B.wtn_rst_guni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.xtr_rst_guni_aom_EPD_2 <- dplyr::select(B.xtr_rst_guni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.btw_rst_uuni_aom_EPD_2 <- dplyr::select(B.btw_rst_uuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.wtn_rst_uuni_aom_EPD_2 <- dplyr::select(B.wtn_rst_uuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.xtr_rst_uuni_aom_EPD_2 <- dplyr::select(B.xtr_rst_uuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.btw_rst_wuni_aom_EPD_2 <- dplyr::select(B.btw_rst_wuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.wtn_rst_wuni_aom_EPD_2 <- dplyr::select(B.wtn_rst_wuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.xtr_rst_wuni_aom_EPD_2 <- dplyr::select(B.xtr_rst_wuni_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.btw_rst_atch_aom_EPD_2 <- dplyr::select(B.btw_rst_atch_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.wtn_rst_atch_aom_EPD_2 <- dplyr::select(B.wtn_rst_atch_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())
B.xtr_rst_atch_aom_EPD_2 <- dplyr::select(B.xtr_rst_atch_aom_EPD_1, 
                                          metric, category, Comparison, 
                                          dplyr::everything())

B.btw_rst_guni_aom_EPD <- B.btw_rst_guni_aom_EPD_2[
  match(B.btw_comp, B.btw_rst_guni_aom_EPD_2$Comparison), ]
B.wtn_rst_guni_aom_EPD <- B.wtn_rst_guni_aom_EPD_2[
  match(B.wtn_comp, B.wtn_rst_guni_aom_EPD_2$Comparison), ]
B.btw_rst_uuni_aom_EPD <- B.btw_rst_uuni_aom_EPD_2[
  match(B.btw_comp, B.btw_rst_uuni_aom_EPD_2$Comparison), ]
B.wtn_rst_uuni_aom_EPD <- B.wtn_rst_uuni_aom_EPD_2[
  match(B.wtn_comp, B.wtn_rst_uuni_aom_EPD_2$Comparison), ]
B.btw_rst_wuni_aom_EPD <- B.btw_rst_wuni_aom_EPD_2[
  match(B.btw_comp, B.btw_rst_wuni_aom_EPD_2$Comparison), ]
B.wtn_rst_wuni_aom_EPD <- B.wtn_rst_wuni_aom_EPD_2[
  match(B.wtn_comp, B.wtn_rst_wuni_aom_EPD_2$Comparison), ]
B.btw_rst_atch_aom_EPD <- B.btw_rst_atch_aom_EPD_2[
  match(B.btw_comp, B.btw_rst_atch_aom_EPD_2$Comparison), ]
B.wtn_rst_atch_aom_EPD <- B.wtn_rst_atch_aom_EPD_2[
  match(B.wtn_comp, B.wtn_rst_atch_aom_EPD_2$Comparison), ]
B.xtr_rst_guni_aom_EPD <- B.xtr_rst_guni_aom_EPD_2 # rename for consistent _#
B.xtr_rst_uuni_aom_EPD <- B.xtr_rst_uuni_aom_EPD_2 # rename for consistent _#
B.xtr_rst_wuni_aom_EPD <- B.xtr_rst_wuni_aom_EPD_2 # rename for consistent _#
B.xtr_rst_atch_aom_EPD <- B.xtr_rst_atch_aom_EPD_2 # rename for consistent _#
B.rst_guni_aom_EPD <- rbind(B.btw_rst_guni_aom_EPD, B.wtn_rst_guni_aom_EPD, 
                            B.xtr_rst_guni_aom_EPD)
B.rst_uuni_aom_EPD <- rbind(B.btw_rst_uuni_aom_EPD, B.wtn_rst_uuni_aom_EPD, 
                            B.xtr_rst_uuni_aom_EPD)
B.rst_wuni_aom_EPD <- rbind(B.btw_rst_wuni_aom_EPD, B.wtn_rst_wuni_aom_EPD, 
                            B.xtr_rst_wuni_aom_EPD)
B.rst_atch_aom_EPD <- rbind(B.btw_rst_atch_aom_EPD, B.wtn_rst_atch_aom_EPD, 
                            B.xtr_rst_atch_aom_EPD)

B.rst_bdiv_aom_EPD_0 <- rbind(B.rst_guni_aom_EPD, B.rst_uuni_aom_EPD, 
                              B.rst_wuni_aom_EPD, B.rst_atch_aom_EPD)

# for _cln, format the data as follows:
# (1) combine the results together
# (2:3) copy and add column category

B.rst_bdiv_aom_cln_0 <- rbind(B.prm_guni_aom_cln, B.prm_uuni_aom_cln,
                              B.prm_wuni_aom_cln, B.prm_atch_aom_cln)
B.rst_bdiv_aom_cln_1 <- B.rst_bdiv_aom_cln_0
B.rst_bdiv_aom_cln_1$category <- "between"

# for _EPD and _cln, format the data as follows:
# (1) rename columns 'p_raw' and 'p_BH'
# (2) create cols 'sig_' with symbols for significance of P and BH_P values
# (3) coerce cols with symbols for significance of P & BH_P values to character
# (4:5) copy and add new column specifying the statistical test
# (6) reorder columns

B.rst_bdiv_aom_EPD_1 <- dplyr::rename(B.rst_bdiv_aom_EPD_0, P = p_raw,
                                      BH_P = p_BH)
B.rst_bdiv_aom_EPD_2 <- dplyr::mutate(
  B.rst_bdiv_aom_EPD_1, 
  sig_P = symnum(P, cutpoints = cutpts, symbols = symbls, corr = F),
  sig_BH_P = symnum(BH_P, cutpoints = cutpts, symbols = symbls, corr = F))
B.rst_bdiv_aom_EPD_3 <- dplyr::mutate(B.rst_bdiv_aom_EPD_2,
                                      sig_P = as.character(sig_P),
                                      sig_BH_P = as.character(sig_BH_P))
B.rst_bdiv_aom_EPD_4 <- B.rst_bdiv_aom_EPD_3
B.rst_bdiv_aom_EPD_4$Test <- "PERMANOVA (adonis)"
B.rst_bdiv_aom_EPD <- dplyr::select(B.rst_bdiv_aom_EPD_4, metric, category, 
                                    Comparison, Test, R2, P, sig_P, BH_P, 
                                    sig_BH_P) # this object is output

B.rst_bdiv_aom_cln_2 <- dplyr::rename(B.rst_bdiv_aom_cln_1, P = p_raw,
                                      BH_P = p_BH)
B.rst_bdiv_aom_cln_3 <- dplyr::mutate(
  B.rst_bdiv_aom_cln_2, 
  sig_P = symnum(P, cutpoints = cutpts, symbols = symbls, corr = F),
  sig_BH_P = symnum(BH_P, cutpoints = cutpts, symbols = symbls, corr = F))
B.rst_bdiv_aom_cln_4 <- dplyr::mutate(B.rst_bdiv_aom_cln_3,
                                      sig_P = as.character(sig_P),
                                      sig_BH_P = as.character(sig_BH_P))
B.rst_bdiv_aom_cln_5 <- B.rst_bdiv_aom_cln_4
B.rst_bdiv_aom_cln_5$Test <- "PERMANOVA (adonis)"
B.rst_bdiv_aom_cln <- dplyr::select(B.rst_bdiv_aom_cln_5, metric, category, 
                                    Comparison, Test, R2, P, sig_P, BH_P, 
                                    sig_BH_P)

## filter results
# for _EPD, this process occurs as follows:
# (1) format to retain only the btw and wtn comparisons
# (2:3) copy & and coerce column Comparison into a character
# (4:6) copy & add full descriptions for info Diet.Type info in col Comparison

B.rfl_bdiv_aom_EPD_0 <- dplyr::filter(B.rst_bdiv_aom_EPD, 
                                      category == "between" |
                                        category == "within")
B.rfl_bdiv_aom_EPD_1 <- B.rfl_bdiv_aom_EPD_0
B.rfl_bdiv_aom_EPD_1$Comparison <- as.character(B.rfl_bdiv_aom_EPD_1$Comparison)
B.rfl_bdiv_aom_EPD <- B.rfl_bdiv_aom_EPD_1
B.rfl_bdiv_aom_EPD$CompKey <- B.rfl_bdiv_aom_EPD$Comparison
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = ".E", 
                                   replacement = ".cecum",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = ".P", 
                                   replacement = ".proximal colon",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = ".D", 
                                   replacement = ".distal colon",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = "C.", 
                                   replacement = "Control ",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = "R.", 
                                   replacement = "Rice bran ",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)
B.rfl_bdiv_aom_EPD$CompKey <- gsub(pattern = "F.", 
                                   replacement = "Fermented rice bran ",
                                   x = B.rfl_bdiv_aom_EPD$CompKey)

# for _cln, this process occurs as follows:
# (1:2) copy & and coerce column Comparison into a character
# (3:5) copy & add full descriptions for info DietGroupLab info in col Comparison

B.rfl_bdiv_aom_cln_0 <- B.rst_bdiv_aom_cln
B.rfl_bdiv_aom_cln_0$Comparison <- as.character(B.rfl_bdiv_aom_cln_0$Comparison)

B.rfl_bdiv_aom_cln <- B.rfl_bdiv_aom_cln_0
B.rfl_bdiv_aom_cln$CompKey <- B.rfl_bdiv_aom_cln$Comparison
B.rfl_bdiv_aom_cln$CompKey <- gsub(
  pattern = "C vs F", 
  replacement = "Control versus Fermented rice bran (colon combined)",
  x = B.rfl_bdiv_aom_cln$CompKey)
B.rfl_bdiv_aom_cln$CompKey <- gsub(
  pattern = "C vs R", 
  replacement = "Rice bran versus Control (colon combined)",
  x = B.rfl_bdiv_aom_cln$CompKey)
B.rfl_bdiv_aom_cln$CompKey <- gsub(
  pattern = "R vs F", 
  replacement = "Rice bran versus Fermented rice bran (colon combined)",
  x = B.rfl_bdiv_aom_cln$CompKey)

# combine the _rfl data.frames and reorder columns
B.rst_flt_bdiv_aom_0 <- rbind(B.rfl_bdiv_aom_cln, B.rfl_bdiv_aom_EPD)
B.rst_flt_bdiv_aom <- dplyr::select(B.rst_flt_bdiv_aom_0,
                                    metric, category, Comparison, Test,
                                    R2, P, sig_P, BH_P, sig_BH_P, CompKey)

### ************************************
### B - STEP  5 - define plot parameters/customize plot aesthetics ----
### ************************************

# NOTE: for DietGroupLab aesthetics, order = C, R, F
# NOTE: for Diet.Type aesthetics, order = C.P, C.D, R.P, R.D, F.P, F.D
# these orders are defined at the end of "B - STEP 2" and "B - STEP 3"

## shapes, colors, labels
# (1) vector for cecum DietGroupLab shapes
# (2) vector for cecum DietGroupLab colors
# (3) vector for cecum DietGroupLab fill
# (4) vector for colon Diet.Type shapes
# (5) vector for colon Diet.Type colors
# (6) vector for colon Diet.Type fill

B.shp_diet_cec <- c(shp_C, shp_R, shp_F)
B.hex_diet_cec <- c(hex_C[1], hex_R[1], hex_F[1])
B.fil_diet_cec <- c(hex_C[2], hex_R[2], hex_F[2])
B.shp_diet_cln <- c(shp_C, shp_C, shp_R, shp_R, shp_F, shp_F)
B.hex_diet_cln <- c(hex_C[1], hex_C[2], hex_R[1], hex_R[2], hex_F[1], hex_F[2])
B.fil_diet_cln <- c(hex_C[2], hex_C[4], hex_R[2], hex_R[4], hex_F[2], hex_F[4])

## sizing parameters
B.sze_int_lne <- 0.31 # size for origin intercept lines
B.sze_int_crc <- 1.81 # size for origin circle
B.sze_int_sqr <- 0.42 # size for origin square
B.sze_smp_pts <- 1.51 # size for sample points
B.sze_ptl <- 16 # size for plot title (i.e. Type)
B.sze_stl <- 14 # size for plot title (i.e. metric)
B.sze_atl <- 14 # size for axis titles (y)
B.sze_ltx <- 14 # size for legend text

# plot and legend labeling
B.ttl_cec <- text_grob("Cecum", face = "bold", color = greydient[1], 
                       size = B.sze_ptl)
B.ttl_cln <- text_grob("Colon (proximal and distal)", face = "bold",
                       color = greydient[1], size = B.sze_ptl)
B.ttl_guni <- "GUniFrac"
B.ttl_atch <- "Aitchison Distance"

## scalers
B.scl_lgn <- 2.00 # value to scale point size in legend

# custom ggplot theme parameters for PCoA and PCA
B.JacksonP <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  line = element_line(lineend = "square", color = greydient[1]),
  plot.title = element_text(size = B.sze_stl, hjust = 0.5),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(size = B.sze_atl),
  legend.position = "none",
  panel.ontop = T,
  panel.background = element_rect(fill = NA))

# custom ggplot theme parameters for legends
B.PurpleRain <- theme(
  text = element_text(family = fnt, color = greydient[1]),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(fill = NA, color = NA),
  legend.box.margin = margin(1, 1, 2, 1, "mm"),
  legend.key.size = unit(0, "mm"),
  legend.text = element_text(
    margin = margin(1, 2, 1, -1, "mm"), size = B.sze_ltx, 
    color = greydient[1], hjust = 0))

### ************************************
### B - STEP 6a - plotting: create PCoA and PCA frameworks (i.e. skeletons) ----
### ************************************

# create empty plots with origin intercept lines (ensures proper layering)

# guni
B.gpr_seg_guni_aom_E <- ggscatter(data = B.smp_cord_guni_aom_E,
                                  x = "xvar", y = "yvar", shape = 32,
                                  font.family = fnt, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8],
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7],
             size = B.sze_int_sqr)

B.gpr_seg_guni_aom_PD <- ggscatter(data = B.smp_cord_guni_aom_PD,
                                   x = "xvar", y = "yvar", shape = 32,
                                   font.family = fnt, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8],
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7],
             size = B.sze_int_sqr)

# atch
B.gpr_seg_atch_aom_E <- ggscatter(data = B.smp_cord_atch_aom_E,
                                  x = "xvar", y = "yvar", shape = 32,
                                  font.family = fnt, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8],
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7],
             size = B.sze_int_sqr)

B.gpr_seg_atch_aom_PD <- ggscatter(data = B.smp_cord_atch_aom_PD,
                                   x = "xvar", y = "yvar", shape = 32,
                                   font.family = fnt, color = greydient[8]) +
  geom_hline(yintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_vline(xintercept = 0, linetype = "dotted", color = greydient[7],
             size = B.sze_int_lne) +
  geom_point(aes(x = 0, y = 0), shape = 16, stroke = 0, color = greydient[8],
             size = B.sze_int_crc) +
  geom_point(aes(x = 0, y = 0), shape = 15, stroke = 0, color = greydient[7],
             size = B.sze_int_sqr)

### ************************************
### B - STEP 6b - plotting: create PCoA and PCA scatterplots ----
### ************************************

# provide provenance for information gathering at end of STEP 1i:
B.prov_secstep_BS6b <- "Section B - STEP 6b"
B.prov_heading_BS6b <- "plotting: create PCoA and PCA scatterplots"
# NOTE: each of the plots below are output to the vault

# begin plotting
B.gpr_guni_cec <- ggscatter(data = B.smp_cord_guni_aom_E,
                            ggp = B.gpr_seg_guni_aom_E,
                            xlab = B.pc1_guni_aom_E,
                            ylab = B.pc2_guni_aom_E,
                            title = B.ttl_guni,
                            x = "xvar", y = "yvar",
                            #shape = shp_cec,
                            shape = "DietGroupLab",
                            fill = "DietGroupLab",
                            color = "DietGroupLab",
                            size = B.sze_smp_pts,
                            font.family = fnt) +
  scale_shape_manual(values = B.shp_diet_cec) +
  scale_color_manual(values = B.hex_diet_cec) +
  scale_fill_manual(values = B.fil_diet_cec) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_atch_cec <- ggscatter(data = B.smp_cord_atch_aom_E,
                            ggp = B.gpr_seg_atch_aom_E,
                            xlab = B.pc1_atch_aom_E,
                            ylab = B.pc2_atch_aom_E,
                            title = B.ttl_atch,
                            x = "xvar", y = "yvar",
                            shape = "DietGroupLab",
                            fill = "DietGroupLab",
                            color = "DietGroupLab",
                            size = B.sze_smp_pts,
                            font.family = fnt) +
  scale_shape_manual(values = B.shp_diet_cec) +
  scale_color_manual(values = B.hex_diet_cec) +
  scale_fill_manual(values = B.fil_diet_cec) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_guni_cln <- ggscatter(data = B.smp_cord_guni_aom_PD,
                            ggp = B.gpr_seg_guni_aom_PD,
                            xlab = B.pc1_guni_aom_PD,
                            ylab = B.pc2_guni_aom_PD,
                            title = B.ttl_guni,
                            x = "xvar", y = "yvar",
                            shape = "Diet.Type",
                            fill = "Diet.Type",
                            color = "Diet.Type",
                            size = B.sze_smp_pts,
                            font.family = fnt) +
  scale_shape_manual(values = B.shp_diet_cln) +
  scale_color_manual(values = B.hex_diet_cln) +
  scale_fill_manual(values = B.fil_diet_cln) +
  border(color = greydient[1]) +
  B.JacksonP

B.gpr_atch_cln <- ggscatter(data = B.smp_cord_atch_aom_PD,
                            ggp = B.gpr_seg_atch_aom_PD,
                            xlab = B.pc1_atch_aom_PD,
                            ylab = B.pc2_atch_aom_PD,
                            title = B.ttl_atch,
                            x = "xvar", y = "yvar",
                            shape = "Diet.Type",
                            fill = "Diet.Type",
                            color = "Diet.Type",
                            size = B.sze_smp_pts,
                            font.family = fnt) +
  scale_shape_manual(values = B.shp_diet_cln) +
  scale_color_manual(values = B.hex_diet_cln) +
  scale_fill_manual(values = B.fil_diet_cln) +
  border(color = greydient[1]) +
  B.JacksonP

### ************************************
### B - STEP 6c - plotting: create custom legends ----
### ************************************

# DietGroup legend for cecum
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels
B.lgn_cec <- data.frame(group = c("A", "B", "C"), x = 0, y = 0,
                        stringsAsFactors = F)
B.lgn_cec_lab <- c(expression(~Con^spf), 
                   expression(~RB^spf), 
                   expression(~FRB^spf))
# plot
B.gpr_lgn_cec <- ggscatter(data = B.lgn_cec, x = "x", y = "y", 
                           shape = "group", fill = "group", color = "group",
                           size = B.sze_smp_pts * B.scl_lgn, 
                           font.family = fnt) +
  scale_shape_manual(values = B.shp_diet_cec, labels = B.lgn_cec_lab,
                     name = NULL) +
  scale_color_manual(values = B.hex_diet_cec, labels = B.lgn_cec_lab,
                     name = NULL) +
  scale_fill_manual(values = B.fil_diet_cec, labels = B.lgn_cec_lab, 
                    name = NULL) +
  B.PurpleRain

# extract legend and convert from class 'gtable' to 'ggplot'
B.gtb_lgn_cec <- get_legend(B.gpr_lgn_cec)
B.ggp_lgn_cec <- as_ggplot(B.gtb_lgn_cec)

# DietGroup and Type legend for colon
# create data.frames with arbitrary x and y values that will be plotted ...
# ... and define legend breaks and labels
B.lgn_cln <- data.frame(group = c("A", "B", "C", "D", "E", "F"), x = 0, y = 0,
                        stringsAsFactors = F)
B.lgn_cln_lab <- c(expression(~Con^spf), 
                   expression(~RB^spf), 
                   expression(~FRB^spf),
                   " ", "proximal", "distal")

B.lgn_shp_cln <- c(B.shp_diet_cec, 32, 22, 22)
B.lgn_hex_cln <- c(B.fil_diet_cec, greydient[8], greydient[1], greydient[1])
B.lgn_fil_cln <- c(B.fil_diet_cec, greydient[8], greydient[3], greydient[7])

# plot
B.gpr_lgn_cln <- ggscatter(data = B.lgn_cln, x = "x", y = "y", 
                           shape = "group", fill = "group", color = "group",
                           size = B.sze_smp_pts * B.scl_lgn, 
                           font.family = fnt) +
  scale_shape_manual(values = B.lgn_shp_cln, labels = B.lgn_cln_lab, 
                     name = NULL) +
  scale_color_manual(values = B.lgn_hex_cln, labels = B.lgn_cln_lab, 
                     name = NULL) +
  scale_fill_manual(values = B.lgn_fil_cln, labels = B.lgn_cln_lab, 
                    name = NULL) +
  guides(shape = guide_legend(nrow = 1)) +
  B.PurpleRain

# extract legend and convert from class 'gtable' to 'ggplot'
B.gtb_lgn_cln <- get_legend(B.gpr_lgn_cln)
B.ggp_lgn_cln <- as_ggplot(B.gtb_lgn_cln)

### ************************************
### B - STEP 6d - plotting: arrange plots and format ----
### ************************************

# provide provenance for information gathering at end of STEP 1i:
B.prov_secstep_BS6d <- "Section B - STEP 6d"
B.prov_heading_BS6d <- "plotting: arrange plots and format"
# NOTE: each of the plots below are output to the vault

# the arrangement process occurs as follows:
# (1) arrange by Type
# (2) arrange with legend
# (3) define Type labels
# (4) annotate arrangements with Type labels and panel labels
# rinse and repeat

B.gga_cec_0 <- ggarrange(B.gpr_guni_cec, B.gpr_atch_cec, labels = NULL,
                         ncol = 2, nrow = 1, align = "hv")
B.gga_cec_1 <- ggarrange(B.gga_cec_0, B.ggp_lgn_cec, labels = NULL,
                         ncol = 1, nrow = 2, heights = c(4, 0.5))
B.gga_cec <- annotate_figure(p = B.gga_cec_1, top = B.ttl_cec,
                             fig.lab = "C", fig.lab.face = "bold",
                             fig.lab.size = 30)

B.gga_cln_0 <- ggarrange(B.gpr_guni_cln, B.gpr_atch_cln, labels = NULL,
                         ncol = 2, nrow = 1, align = "hv")
B.gga_cln_1 <- ggarrange(B.gga_cln_0, B.ggp_lgn_cln, labels = NULL,
                         ncol = 1, nrow = 2, heights = c(4, 0.5))
B.gga_cln <- annotate_figure(p = B.gga_cln_1, top = B.ttl_cln,
                             fig.lab = "D", fig.lab.face = "bold",
                             fig.lab.size = 30)

### ************************************
### B - WRITE OUTPUTS ----
### ************************************

# provenance for outputs to the vault
B.prov_output1 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_rst_bdiv_EPD,
                             "object" = B.prov_outobj1_BS4b,
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS4b,
                             "heading" = B.prov_heading_BS4b,
                             stringsAsFactors = F)
B.prov_output2 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_rst_flt_bdiv,
                             "object" = B.prov_outobj2_BS4b,
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS4b,
                             "heading" = B.prov_heading_BS4b,
                             stringsAsFactors = F)
B.prov_output3 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gpr_guni_cec,
                             "object" = "B.gpr_guni_cec",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6b,
                             "heading" = B.prov_heading_BS6b,
                             stringsAsFactors = F)
B.prov_output4 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gpr_atch_cec,
                             "object" = "B.gpr_atch_cec",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6b,
                             "heading" = B.prov_heading_BS6b,
                             stringsAsFactors = F)
B.prov_output5 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gpr_guni_cln,
                             "object" = "B.gpr_guni_cln",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6b,
                             "heading" = B.prov_heading_BS6b,
                             stringsAsFactors = F)
B.prov_output6 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gpr_atch_cln,
                             "object" = "B.gpr_atch_cln",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6b,
                             "heading" = B.prov_heading_BS6b,
                             stringsAsFactors = F)
B.prov_output7 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gga_cec,
                             "object" = "B.gga_cec",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6d,
                             "heading" = B.prov_secstep_BS6d,
                             stringsAsFactors = F)
B.prov_output8 <- data.frame("info" = "provenance for output",
                             "path" = B.ofv_gga_cln,
                             "object" = "B.gga_cln",
                             "script" = name_scrpt,
                             "section" = B.prov_secstep_BS6d,
                             "heading" = B.prov_secstep_BS6d,
                             stringsAsFactors = F)
B.prov <- rbind(B.prov_output1, B.prov_output2, B.prov_output3, B.prov_output4,
                B.prov_output5, B.prov_output6, B.prov_output7, B.prov_output8)

# outputs to the vault
write.table(sep = "\t", row.names = F, x = B.rst_bdiv_aom_EPD, 
            file = B.ofv_rst_bdiv_EPD)
write.table(sep = "\t", row.names = F, x = B.rst_flt_bdiv_aom, 
            file = B.ofv_rst_flt_bdiv)

write.table(sep = "\t", row.names = F, x = B.prov, file = B.ofv_prov)

ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 75,
       filename = B.ofv_gpr_guni_cec, plot = B.gpr_guni_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 75,
       filename = B.ofv_gpr_atch_cec, plot = B.gpr_atch_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 75,
       filename = B.ofv_gpr_guni_cln, plot = B.gpr_guni_cln)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 75,
       filename = B.ofv_gpr_atch_cln, plot = B.gpr_atch_cln)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = B.ofv_gga_cec, plot = B.gga_cec)
ggsave(device = "pdf", dpi = 600, units = "mm", width = wid_dbl, height = 90,
       filename = B.ofv_gga_cln, plot = B.gga_cln)

B.obj <- ls(pattern = "B.")
B.lst <- c(B.obj[grep(pattern = "B.", x = B.obj, ignore.case = F, fixed = T)],
           PREFACE.lst, COSMOS, B.obj_from_ST)
save(list = B.lst, file = B.ofv_wksp)
