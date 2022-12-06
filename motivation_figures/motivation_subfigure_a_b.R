library(ggplot2)
library(GenomicRanges)
library(dplyr)

dt <- read.table("segdups.tsv",  header=TRUE)

dt$size <- dt$chromEnd - dt$chromStart

ggplot(dt, aes(x = fracMatch, y = size)) + 
  geom_density_2d_filled(show.legend = FALSE) + 
  scale_fill_brewer() +
  theme_minimal() +
  theme(legend.text = element_text(size = 14),   
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        strip.text.y = element_blank()) +
  scale_y_log10(  
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),  limits=c(10^3, 1.8*10^6)) + 
  annotation_logticks(sides = "l") +
  xlab("Size of SD (log scale)") + 
  ylab("Percent identity") 
  
get_stats <- function(dt, inside) {

  dt <- dt[(dt$peri==0 & dt$telo==0 &dt$acro==0) == inside, ]
  
  gr <- GRanges(
    seqnames=dt$chrom,
    ranges=IRanges(start = dt$chromStart, end = dt$chromEnd),
  )
  
  gr_d <- disjoin(gr)
  overlaps_no <- countOverlaps(gr_d, gr)
  gr_width <- width(gr_d)
  
  overlap_data <- data.frame(overlaps_no = overlaps_no, width=gr_width)
  results <- overlap_data  %>% group_by(overlaps_no) %>%  summarize(stat=sum( width))
  results$inside <- inside
  return (results)
}

dt <- dt[dt$size > 5000,]
stats_inside <- get_stats(dt, TRUE)
stats_outside <- get_stats(dt, FALSE)

stats <- rbind(stats_inside, stats_outside)

ggplot(stats, aes(x = overlaps_no, y = stat, fill = inside)) +
  geom_bar(stat = "identity") +
  xlab("Number of overlapped SDs") + ylab("Total length of SDs") +
  scale_y_log10(  
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  labs(fill="Genomic Region") +
  scale_fill_manual(values = c("FALSE" = "#f7a716",
                               "TRUE" = "#31479e"),
                    labels=c( "non-interstitial", "interstitial")) +
  annotation_logticks(sides = "l") +
  theme_minimal() +
  theme(legend.position="bottom", 
        legend.text = element_text(size = 14),   
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        strip.text.y = element_blank()) +
  guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(inside~.) +
  xlim(0,150) 
  
