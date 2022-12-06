library(ggplot2)
library(GenomicRanges)
library(dplyr)

dt <- read.table("segdups.tsv",  header=TRUE)
dt$size <- dt$chromEnd - dt$chromStart
dt <- dt[dt$size > 5000,]


get_stats <- function(dt, treshold) {
  
  dt <- dt[ dt$fracMatch >= treshold & (dt$peri==0 & dt$telo==0 &dt$acro==0), ]
  
  gr <- GRanges(
    seqnames=dt$chrom,
    ranges=IRanges(start = dt$chromStart, end = dt$chromEnd)
  )
  
  gr_d <- disjoin(gr)
  overlaps_no <- countOverlaps(gr_d, gr)
  gr_width <- width(gr_d)
  
  zeros_data <- data.frame(overlaps_no = 1:500, width=0)
  overlap_data <- data.frame(overlaps_no = overlaps_no, width=gr_width)
  overlap_data <- rbind(zeros_data,  overlap_data )
  
  results <- overlap_data  %>% group_by(overlaps_no) %>%  summarize(stat_sum=sum(width))
  results <- results %>%arrange(desc(overlaps_no)) %>% mutate(stat=cumsum(stat_sum))

  results$treshold <- treshold
  return (results)
}

v <- c(0.9, 0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 1)

d_all <- data_frame()

for (i in 1:(length(v)-1)) {
  d <- get_stats(dt, v[i])
  d_all <- rbind(d, d_all)
} 

dark_blue = "#31479e"
middle_blue = "#4d6ab5"
light_blue ="#6e91c8"
green = "#b8dcb0"
yellow = "#efed1d"
orange = "#f7a716"
claret = "#94171b"
red = "#e62427"

colors = rev(c(claret,red,   orange,yellow, light_blue, middle_blue, dark_blue))

ggplot(d_all, aes(x = overlaps_no, y = stat,  colour=factor(treshold)))  + 
  geom_ribbon(data = d_all, aes(x= overlaps_no, ymin =  0, ymax =  stat, fill = factor(treshold))) +
  scale_y_log10(  
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)), 
    limits=c(10^5, 1.1 * 10^8)) + 
  annotation_logticks(sides = "l") +
  xlim(c(0, 50)) +
  theme_minimal() +
  # labs(colour="Identity treshold") +  
  # xlab("Number of overlapped SDs") + 
  # ylab("Total length of SDs (log scale)") + 
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom", 
        legend.text = element_text(size = 14),   
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 16), 
        axis.text.y = element_text(size = 14)) 
