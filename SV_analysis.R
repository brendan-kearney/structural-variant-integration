#install.packages("devtools")
#devtools::install_github("venyao/intansv")
#install.packages("openxlsx")
#install.packages("ggplot2")

library(intansv)
library(openxlsx)
library(ggplot2)

# User Configurable Variables
file_directory <- "data/example_data"
annotation_directory <- "data/annotation"
output_directory <- "outputs"

# File paths
delly_file_name <- "acgc170.delly.vcf"
breakdancer_file_name <- "acgc170.breakdancer.filtered.out"
lumpy_file_name <- "acgc170.lumpy.vcf"
pindel_dir_name <- "pindel_out"
annotation_file <- file.path(annotation_directory, "hg38_annotation.txt")
genome_file <- file.path(annotation_directory, "hg38_sizes.copy.txt")

# Derived file paths
delly.file.path <- file.path(file_directory, delly_file_name)
breakdancer.file.path <- file.path(file_directory, breakdancer_file_name)
lumpy.file.path <- file.path(file_directory, lumpy_file_name)
pindel.dir.path <- file.path(file_directory, pindel_dir_name)

# DELLY is able to predict deletions, inversions and duplications
delly <- readDelly(delly.file.path)
View(delly)

# Pindel is able to predict deletions, inversions and duplications.
pindel <- readPindel(pindel.dir.path)
View(pindel)

# BreakDancer is able to predict deletions and inversions. 

#breakdancer_data <- read.table(breakdancer.file.path, comment.char="#", header=FALSE,sep="\t", stringsAsFactors=FALSE)
breakdancer <- readBreakDancer(breakdancer.file.path)
View(breakdancer)

# CNVnator/CNVpytor is able to predict deletions and duplications. 
# To do: match output format with what's required by intansv
#cnvnator <- readCnvnator(cnvnator.dir.path)

# SVseq2 is able to predict deletions.
# To do: finish chr2 run, import

#svseq <- readSvseq(svseq.dir.path)

lumpy <- readLumpy(lumpy.file.path)
str(lumpy)
sv_all_methods <- methodsMerge(breakdancer,delly,pindel,lumpy)

anno.file.path <- file.path(file_directory, "hg38_annotation.txt")
hg38_gff <- read.table(anno.file.path, head=TRUE)
sv_all_methods.anno <- llply(sv_all_methods,svAnnotation,genomeAnnotation=hg38_gff)
View(sv_all_methods.anno)

genome.file.path <- file.path("/Users/btk20/Documents/SV_analysis/hg38_sizes.copy.txt")
genome <- read.table(genome.file.path, head=TRUE, as.is=TRUE)
plotRegion(sv_all_methods,hg38_gff, "chr1",1,20000000)
write.xlsx(sv_all_methods.anno, file="/Users/btk20/Documents/SV_analysis/sv_all_methods.anno.xlsx")

# Count the number of methods in each data frame (del, dup, inv)
sv_all_methods$del$num_methods <- sapply(strsplit(sv_all_methods$del$methods, ":"), length)
sv_all_methods$dup$num_methods <- sapply(strsplit(sv_all_methods$dup$methods, ":"), length)
sv_all_methods$inv$num_methods <- sapply(strsplit(sv_all_methods$inv$methods, ":"), length)

# Combine the data into one data frame with SV_Type and num_methods
sv_data <- rbind(
  data.frame(SV_Type = "del", num_methods = sv_all_methods$del$num_methods),
  data.frame(SV_Type = "dup", num_methods = sv_all_methods$dup$num_methods),
  data.frame(SV_Type = "inv", num_methods = sv_all_methods$inv$num_methods)
)

# Create a summary data frame that counts the number of occurrences by SV_Type and num_methods
sv_summary <- as.data.frame(table(sv_data$SV_Type, sv_data$num_methods))
names(sv_summary) <- c("SV_Type", "num_methods", "Count")
sv_summary$num_methods <- as.factor(sv_summary$num_methods)  # Ensure num_methods is a factor
#295, 305, 306, 323, 446, 455, 470

# Create the stacked bar chart
ggplot(sv_summary, aes(x = SV_Type, y = Count, fill = num_methods)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Observations for Each SV Type by Number of Methods",
       x = "Structural Variation Type",
       y = "Number of Observations",
       fill = "Number of Methods") +
  scale_fill_brewer(palette = "Set1")
