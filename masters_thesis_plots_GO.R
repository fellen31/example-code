library(tidyverse)
library(data.table)
library(ggpubr)
library(GOfuncR)
library(MetBrewer)
library(png)

################
### Figure 2 ###
################

calanus_img <- readPNG("~/calanus2.png")
sampling_img <- readPNG("~/sampling.png")

fig2A <- ggplot() + 
  background_image(calanus_img) +
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

fig2B <- ggplot() + 
  background_image(sampling_img) +
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

# Set path to GMAP psl files
path = "/vaults/nvme2_2tb/felix/projects/tmp"

# Set filenames 
files = dir(path, pattern = "*.psl", recursive = TRUE)
# Read files into df

# Read all files 
psl <- tibble(file = files) %>%
  mutate(file_contents = map(file, 
                             ~ read_delim(file.path(path, .), 
                                          delim = '\t', 
                                          skip = 0, 
                                          col_names = c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"), 
                                          col_types = cols(.default = "?", tStarts = "c")))) %>% unnest(file_contents)

# Filter psl files for plotting
psl_pretty_print <- psl %>%
  # Extract info from psl file names 
  mutate(assembly = gsub("\\.index.*", "", file), 
         sample = gsub(".*index\\.", "", file),
         sample = gsub("\\.trinity.*", "", sample), 
         file = gsub("\\.trinity.*", "", file)) %>%
  # Choose versions of assembly to plot
  filter(assembly %in% c("assembly.fasta", "assembly.medaka.fasta", "flye_andreas.medaka.pilon_no_ncbi.fasta"))

# Change names for prettier plotting
psl_pretty_print$assembly[which(psl_pretty_print$assembly == "assembly.fasta")] <- "Flye (Racon)"
psl_pretty_print$assembly[which(psl_pretty_print$assembly == "assembly.medaka.fasta")] <- "ONT (Medaka)"
psl_pretty_print$assembly[which(psl_pretty_print$assembly == "flye_andreas.medaka.pilon_no_ncbi.fasta")] <- "RNA (Pilon)"

# Save to be able to subset text in plot
inserts_plot_data <- psl_pretty_print %>%
  # Use just one sample to plot
  filter(sample == "UI-3062-CalhypRNA1_S11_L002") %>%
  mutate(identity = matches/(qEnd-qStart), percent_aligned = (qEnd - qStart)/qSize) %>% 
  filter(identity > 0.8,
         percent_aligned > 0.8) %>%
  # Plot levels instead of numbers
  mutate(qNumInsert = factor(qNumInsert, levels = c(3:0)),
         assembly = factor(assembly, levels = c("Flye (Racon)", 
                                                "ONT (Medaka)", 
                                                "RNA (Pilon)"))) %>%
  # To be able to add everything over 3 inserts into level 3+
  mutate(Inserts = plyr::revalue(qNumInsert, c("3"="3+"))) %>%
  # Remove NA values
  filter(!is.na(Inserts))

# Number of inserts plot 
fig2C <- ggplot(inserts_plot_data, aes(x = assembly)) +
  geom_bar(aes(fill = Inserts))  +
  scale_fill_manual(values = met.brewer("Egypt", 5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  ylab("") +
  geom_text(data = subset(inserts_plot_data, Inserts %in% c("0","1")), stat='count', aes(label=scales::comma(..count..), fill = Inserts), vjust=-1 , position = position_stack(
    vjust = .3, 
    reverse = FALSE), size = 3, vjust = -0.5, alpha=0.5)

# Add data for busco plot
busco_plot_data <- tibble(version = c("Flye (Racon)",
                                      "ONT (Medaka)",
                                      "RNA (Pilon)"),
                          busco = c("C:91.6%[S:66.3%,D:25.3%],F:2.9%,M:5.5%,n:1066",
                                    "C:92.0%[S:66.5%,D:25.5%],F:3.1%,M:4.9%,n:1066",
                                    "C:94.3%[S:64.1%,D:30.2%],F:1.3%,M:4.4%,n:1066")) %>%
  # Disentangle BUSCO syntax 
  mutate(Single = as.numeric(gsub(".*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*", "\\1", busco))/100,
         Duplicated = as.numeric(gsub(".*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*", "\\2", busco))/100,
         Fragmented = as.numeric(gsub(".*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*?([0-9]+\\.[0-9]+).*", "\\3", busco))/100
  ) %>%
  select(!c(busco)) %>%
  pivot_longer(-version) %>%
  # Change names for plotting
  mutate(BUSCO = factor(name, levels = c("Fragmented", "Duplicated", "Single")),
         version = factor(version, levels = c("Flye (Racon)",
                                              "ONT (Medaka)",
                                              "RNA (Pilon)"
         )))

# Plot BUSCO
fig2D <- ggplot(busco_plot_data, aes(x = version, y = value, fill = BUSCO)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = (met.brewer("Egypt", 3))) +
  geom_text(data=subset(busco_plot_data, value >0.1),aes(label=scales::percent(value)), position = position_stack(
    vjust = 0.5, 
    reverse = FALSE), size = 3, vjust = 0, alpha=0.5) +
  xlab("") + 
  ylab("") +
  scale_y_continuous(labels = scales::percent)

# Plot Figure 1 
ggarrange(fig1A + labs(tag = "A"),fig1B + labs(tag = "B"),fig1C + labs(tag = "C"), fig1D + labs(tag = "D"))

################
### Figure 3 ###
################

# Read depth data and summarize length data per contig as well
depth <- fread("/vaults/nvme2_2tb/felix/projects/calanus/data/external/assembly/flye/depth_histogram.tsv")
setnames(depth, "V1", "contig")
setkey(depth, contig)
lengths <- depth[,.(length=sum(V3)), by = contig]
setkey(lengths, contig)

# Reads SVs 
combisv <- fread("/vaults/external/seagate_8tb_1/felix/projects/calanus/data/external/assembly/flye/sv/combisv/assembly.fasta.minimap2.sam.no_secondary.sorted.md.min_support_2_min_qual_2.min_coverage_2.vcf")
setnames(combisv, "#CHROM", "contig")
idxstats <- read_tsv("/vaults/nvme2_2tb/felix/projects/calanus/data/external/assembly/flye/assembly.fasta.minimap2.sam.no_secondary.sorted.bam.idxstats.tsv", col_names = c("contig", "length", "mapped_reads", "unmapped_reads"))

# Unnest INFO column 
unnested_combisv <- tibble(combisv) %>% 
  mutate(INFO = str_split(INFO, ";")) %>%
  unnest(cols = INFO) %>%
  mutate(LEFT = gsub("=.*", "", INFO),
         RIGHT = gsub(".*=", "", INFO)) %>%
  select(!INFO) %>%
  pivot_wider(names_from = LEFT, values_from = RIGHT)

combined <- idxstats %>% full_join(unnested_combisv) %>%
  mutate(SVLEN = as.numeric(SVLEN))

# Read repeats from Red
red <- fread("/vaults/nvme2_2tb/felix/projects/calanus/data/interim/contig_ends/repeats/racon/Red-AllRepeats_coordinates_RaconAssembly.txt.bed")
red <- red %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
red <- reduce(red) %>% as_tibble() %>% mutate(middle = end-start)

# Divide mapping depth into short and long contigs
short_contigs <- depth %>%
  group_by(contig) %>%
  summarize(length = sum(V3), sum_depth = sum(V2*V3)) %>%
  mutate(mean_depth = sum_depth / length) %>%
  as_tibble() %>%
  mutate(short = ifelse(length < 12500, "<12,500 bp", ">=12,500 bp"))
short_contigs <- short_contigs %>% mutate(short = factor(short_contigs$short, levels = c("<12,500 bp", ">=12,500 bp")))

# Divide SVs into short and long contigs
svtypes <- c("DEL", "DUP", "INV", "INS")
short_svs <- combined %>% filter(SVTYPE %in% svtypes) %>% group_by(contig) %>% summarize(sum_width = sum(SVLEN), length) %>%
  filter(!is.na(sum_width)) %>% distinct() %>% mutate(percentage = sum_width / length) %>% mutate(short = ifelse(length < 12500, "<12,500 bp", ">=12,500 bp"))
short_svs <- short_svs %>% ungroup() %>% mutate(short = factor(short_svs$short, levels = c("<12,500 bp", ">=12,500 bp")))

# Divide repeats short and long contigs
short_red <- red %>% inner_join(., lengths, by = c(seqnames = "contig")) %>% group_by(seqnames) %>% summarize(sum_width = sum(width), length) %>%
  filter(!is.na(sum_width)) %>% distinct() %>% mutate(percentage = sum_width / length) %>% mutate(short = ifelse(length < 12500, "<12,500 bp", ">=12,500 bp"))
short_red <- short_red %>% ungroup() %>% mutate(short = factor(short_red$short, levels = c("<12,500 bp", ">=12,500 bp")))

fig3A <- gghistogram(lengths, x = "length")  +
  labs(tag = "A") +   scale_x_continuous(labels = scales::comma)

fig3B <- ggboxplot(short_contigs, x = "short", y = "mean_depth") + yscale("log2", .format = TRUE) + stat_compare_means(method = "t.test") +
  labs(tag = "B") +  scale_y_log10(labels = scales::comma) + xlab("") + ylab("mean depth")

fig3C <- ggboxplot(short_svs, x = "short", y = "percentage") + stat_compare_means(method = "t.test") +
  labs(tag = "C")+  scale_y_log10(labels = scales::percent) + xlab("")

fig3D <- ggboxplot(short_red, x = "short", y = "percentage") + stat_compare_means(method = "t.test") +
  labs(tag = "D")+  scale_y_continuous(labels = scales::percent) + xlab("")

# Plot Figure 3 
ggarrange(fig3A,fig3B,fig3C,fig3D)

################
### Figure 4 ###
################

# Repeats:
# Set directory 
path = "/vaults/nvme2_2tb/felix/projects/calanus/data/interim/contig_ends/repeats/racon"
# Set filenames 
files = dir(path, pattern = "*.25kb.tsv")
filesrev = dir(path, pattern = "*.25kb.rev.tsv")
# Read files into df
c <- tibble(file = files) %>%
  mutate(file_contents = map(file, 
                             ~ read_delim(file.path(path, .), 
                                          delim = " ", 
                                          skip = 0, 
                                          col_names = c("POS", "COUNT", "TOTAL", "FREQ")))) %>%
  unnest(file_contents)

crev <- tibble(file = filesrev) %>%
  mutate(file_contents = map(file, 
                             ~ read_delim(file.path(path, .), 
                                          delim = " ", 
                                          skip = 0, 
                                          col_names = c("POS", "COUNT", "TOTAL", "FREQ")))) %>%
  unnest(file_contents)

# Change names 
c$file[which(c$file == "Red-AllRepeats_coordinates_RaconAssembly.txt.bed.assembly.fasta.25kb.tsv")] <- "Red"
c$file[which(c$file == "SciRoKo-short-repeat-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.tsv")] <- "Short repeats"
c$file[which(c$file == "SciRoKo-mononucleotide-repeat-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.tsv")] <- "Mononucleotide repeats"
c$file[which(c$file == "TideHunter-LTR-coordinates_RaconAssembly.txt.bed.assembly.fasta.25kb.tsv")] <- "Long Tandem Repeats"
c$file[which(c$file == "TransposonPSI-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.tsv")] <- "Transposon ORFs"

crev$file[which(crev$file == "Red-AllRepeats_coordinates_RaconAssembly.txt.bed.assembly.fasta.25kb.rev.tsv")] <- "Red"
crev$file[which(crev$file == "SciRoKo-short-repeat-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.rev.tsv")] <- "Short repeats"
crev$file[which(crev$file == "SciRoKo-mononucleotide-repeat-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.rev.tsv")] <- "Mononucleotide repeats"
crev$file[which(crev$file == "TideHunter-LTR-coordinates_RaconAssembly.txt.bed.assembly.fasta.25kb.rev.tsv")] <- "Long Tandem Repeats"
crev$file[which(crev$file == "TransposonPSI-coordinates-RaconAssembly.txt.bed.assembly.fasta.25kb.rev.tsv")] <- "Transposon ORFs"

# Plot 1 
fig4A <- tibble(c) %>% full_join(crev, by = c(file = "file", POS = "POS")) %>% mutate(FREQ = (FREQ.x + FREQ.y)/2) %>%
  filter(!file %in% c("Red", "Mononucleotide repeats")) %>%
  mutate(Repeat = gsub("\\.bed.*", "", file)) %>%
  ggplot(aes(POS, FREQ, color = Repeat)) +
  geom_line(size = 1) +
  xlab("POS") +
  scale_x_continuous(labels = scales::comma, limits = c(0,12500)) +
  scale_color_manual(values = met.brewer("Egypt",4)) +
  theme(plot.margin=unit(c(1,0,1,1),"cm")) +
  xlab("Position (bp)") +
  ylab("Per base frequency (%)") +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.12))

#SVs:
# Set directory 
path = "/vaults/nvme2_2tb/felix/projects/calanus/data/interim/contig_ends/svs"
# Set filenames 
files = dir(path, pattern = "*.25kb.tsv")
filesrev = dir(path, pattern = "*.25kb.rev.tsv")
# Read files into df

c <- tibble(file = files) %>%
  mutate(file_contents = map(file, 
                             ~ read_delim(file.path(path, .), 
                                          delim = " ", 
                                          skip = 0, 
                                          col_names = c("POS", "COUNT", "TOTAL", "FREQ")))) %>%
  unnest(file_contents)

crev <- tibble(file = filesrev) %>%
  mutate(file_contents = map(file, 
                             ~ read_delim(file.path(path, .), 
                                          delim = " ", 
                                          skip = 0, 
                                          col_names = c("POS", "COUNT", "TOTAL", "FREQ")))) %>%
  unnest(file_contents)
  
crev %>% select(file) %>% distinct()
c$file[which(c$file == "combisv.vcf.DEL.bed.assembly.fasta.25kb.tsv")] <- "Deletions"
c$file[which(c$file == "combisv.vcf.DUP.bed.assembly.fasta.25kb.tsv")] <- "Duplications"
c$file[which(c$file == "combisv.vcf.INS.bed.assembly.fasta.25kb.tsv")] <- "Insertions"
c$file[which(c$file == "combisv.vcf.INV.bed.assembly.fasta.25kb.tsv")] <- "Inversions"


crev$file[which(crev$file == "combisv.vcf.DEL.bed.assembly.fasta.25kb.rev.tsv")] <- "Deletions"
crev$file[which(crev$file == "combisv.vcf.DUP.bed.assembly.fasta.25kb.rev.tsv")] <- "Duplications"
crev$file[which(crev$file == "combisv.vcf.INS.bed.assembly.fasta.25kb.rev.tsv")] <- "Insertions"
crev$file[which(crev$file == "combisv.vcf.INV.bed.assembly.fasta.25kb.rev.tsv")] <- "Inversions"

fig4B <- tibble(c) %>% full_join(crev, by = c(file = "file", POS = "POS")) %>% mutate(FREQ = (FREQ.x + FREQ.y)/2) %>%
  filter(file %in% c("Deletions", "Duplications", "Insertions", "Inversions")) %>%
  # Shorten filename
  mutate(SV = file) %>%
  ggplot(aes(POS, FREQ, color = SV)) +
  geom_line(size = 1) +
  xlab("POS") +
  scale_x_continuous(labels = scales::comma, limits = c(0,12500)) +
  scale_color_manual(values = met.brewer("Egypt",4)) +
  theme(plot.margin=unit(c(1,0,1,1),"cm")) +
  xlab("Position (bp)") +
  ylab("Per base frequency (%)") +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.02))

# Look at mean depth for contig ends 
depths <- read_delim("~/projects/depth/out", col_names = c("base", "depth", "count"), delim = "\t")
depths_rev <- read_delim("~/projects/depth/out_rev", col_names = c("base", "depth", "count"), delim = "\t")

depths <- depths %>% mutate(mean_depth = depth/count) 
depths_rev <- depths_rev %>% mutate(mean_depth = depth/count) 

fig4C <- depths %>% full_join(depths_rev, by = c(base = "base")) %>% mutate(mean_depth = (mean_depth.x + mean_depth.y)/2) %>% ggplot(aes(x=base, y=mean_depth)) + geom_line() + scale_x_continuous(labels = scales::comma) + xlim(0,12500) + ylab("Mean per base depth") +
  xlab("Position (bp)") +
  ylim(0,20)

# Plot Figure 4 
ggarrange(fig4A + labs(tag = "A"), fig4B + labs(tag = "B"), fig4C + labs(tag = "C"))

###############################
### GO enrichment analysis ####
###############################

# Blast results
mm_dups <- read_tsv("/vaults/nvme2_2tb/felix/projects/calanus/data/interim/blast_duplicate_genes/multi_mapping_genes_best_pairs.combined.bastp", col_names = c("tran", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# Filter SM genes
single_mapping_genes <- psl %>% mutate(identity = matches/(qEnd-qStart), percent_aligned = (qEnd - qStart)/qSize) %>% 
  filter(file == "assembly.fasta.index.UI-3062-CalhypRNA1_S11_L002.trinity.Trinity.fasta.longest_isoforms.fasta.transdecoder.cds.psl") %>%
  filter(identity > 0.8,
         percent_aligned > 0.8) %>%
  group_by(qName) %>% mutate(n = n()) %>% filter(n == 1) 

# Filter MM genes
multi_mapping_genes <- psl %>% mutate(identity = matches/(qEnd-qStart), percent_aligned = (qEnd - qStart)/qSize) %>% 
  filter(file == "assembly.fasta.index.UI-3062-CalhypRNA1_S11_L002.trinity.Trinity.fasta.longest_isoforms.fasta.transdecoder.cds.psl") %>%
  filter(identity > 0.8,
         percent_aligned > 0.8) %>%
  group_by(qName) %>% mutate(n = n()) %>% filter(n > 1)

# Select MM genes, filter out haplotigs 
genes_to_do_go_with <- mm_dups %>% group_by(tran) %>%
  summarize(identity_length = sum(pident/100*(qend-qstart)), sum_qsize = sum(qend-qstart), weighted_ident = identity_length/sum_qsize, qseqid, sseqid) %>% arrange(desc(weighted_ident)) %>% select(weighted_ident) %>% distinct() %>%
  inner_join(., multi_mapping_genes, by = c(tran = "qName")) %>%
  mutate(`Identity` = weighted_ident) %>%
  mutate(`Haplotigs` = ifelse(weighted_ident>0.99, "Potentially", "Probably not")) %>% select(Haplotigs, tran) %>% distinct() %>% filter(Haplotigs == "Probably not") %>% select(tran)

# EnTAP results file header
entap_cols_calanus = c("Query Sequence", "Subject Sequence", "Percent Identical", "Alignment Length", "Mismatches", 
                       "Gap Openings", "Query Start", "Query End", "Subject Start", "Subject End", 
                       "E Value", "Coverage", "Description", "Species", "Taxonomic Lineage", 
                       "Origin Database", "Contaminant", "Informative", "UniProt Database Cross Reference", "UniProt Additional Information", 
                       "UniProt KEGG Terms", "UniProt GO Biological", "UniProt GO Cellular", "UniProt GO Molecular", "EggNOG Seed Ortholog", 
                       "EggNOG Seed E-Value", "EggNOG Seed Score", "EggNOG Predicted Gene", "EggNOG Tax Scope", "EggNOG Tax Scope Max", 
                       "EggNOG Member OGs", "EggNOG Description", "EggNOG KEGG Terms", "EggNOG GO Biological", "EggNOG GO Cellular", 
                       "EggNOG GO Molecular", "EggNOG Protein Domains", "X38")
# Read EnTAP results
entap_calanus <- fread("/vaults/external/seagate_8tb_1/felix/projects/calanus/rna_sequencing/results/entap/UI-3062-CalhypRNA1_S11_L002.trinity.Trinity.fasta.longest_isoforms.fasta.transdecoder.pep/final_results/final_annotations_lvl0.tsv", col.names = entap_cols_calanus, sep = "\t")

# Select the representative transcriptome and combine with entap results 
psl_single_file <- psl %>% filter(file == "assembly.fasta.index.UI-3062-CalhypRNA1_S11_L002.trinity.Trinity.fasta.longest_isoforms.fasta.transdecoder.cds.psl")
full_calanus <- full_join(tibble(psl_single_file), tibble(entap_calanus), by = c(qName = "Query Sequence"))

# Select GO:s from EnTAP
full_go_calanus <- full_calanus %>% 
  # Select GO Columns
  select("qName", 
         "EggNOG GO Biological",
         "EggNOG GO Cellular",
         "EggNOG GO Molecular") %>%
  # Put into single column 
  pivot_longer(!qName, values_to = "GO", names_to = "type") %>%
  # Remove empty
  filter(GO != "") %>%
  filter(!is.na(GO)) %>%
  select(qName, GO) %>%
  # Split into single column for each GO
  mutate(GO = str_split(GO, ",")) %>%
  # Make longer 
  unnest(GO) %>%
  # Filter empty
  filter(GO != "") %>%
  # Remove (L=n)
  mutate(GO = gsub("\\(.*", "", GO)) %>%
  mutate(gene = qName, go_id = GO) %>%
  select(gene, go_id)

# Use all high confindece genes as background
background <- rbind(multi_mapping_genes, single_mapping_genes) %>% select(qName) %>% distinct()
database <- full_go_calanus %>% distinct() %>% filter(gene %in% background$qName)

# Select foreground 
fg <- database %>% 
  filter(gene %in% genes_to_do_go_with$tran) %>%
  mutate(gene = gene, is_candidate = 1) %>% ungroup() %>% select(gene, is_candidate) %>% distinct()

# Do GO analysis 
res <- go_enrich(as.data.frame(fg), annotations = as.data.frame(database))
ref <- refine(res, fwer = 0.05, annotations = data.frame(database))

# GO results
tibble(ref) %>% filter(signif==TRUE) 