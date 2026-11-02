library(bacha)
library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2024)
library(TFBSTools)

# load open chromatin regions
ocr <- read.table("https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed", sep = "\t")
ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")

ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "excluderanges")
excludeGR.hg38 <- query_data[["AH107305"]]

repeat_regions <- ah[['AH99003']]
exclude <- c("Low_complexity", "Satellite", "Simple_repeat")
repeat_regions <- repeat_regions[repeat_regions$repClass %in% exclude]
mask_regions <- c(excludeGR.hg38, repeat_regions)
mask_regions <- reduce(mask_regions)
strand(mask_regions) <- "*"

jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
pwm <- TFBSTools::getMatrixSet(sq24, list(tax_group = "vertebrates", collection = "CORE"))

# subset chr19
ocr <- ocr[seqnames(ocr) == "chr19"]

BachaIndex(
  pwm = pwm[1:5],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  open_regions = ocr,
  mask_regions = mask_regions,
  output_path = "test/hg38_chr19_test.bi",
  window = 50,
  n_background = 10
)
