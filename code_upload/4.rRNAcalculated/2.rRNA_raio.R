#### calculate the rRNA ratio
library(data.table)
library(ggplot2)
library(ggsci)
library(languageserver)
library(dplyr)
library(ggplot2)
library(ggforce)
setwd('/mnt/data3/yiyonghao/NC_paper_rawdata')

# 读取mapping_summarySSU.txt文件
SSU = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_calculate/20k_subset/mapping_summary_SSU.txt')
colnames(SSU) = c('sample','mapped','unmapped')
LSU = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_calculate/20k_subset/mapping_summary_LSU.txt')
colnames(LSU) = c('sample','mapped','unmapped')

total_mapped   <- sum(SSU$mapped) + sum(LSU$mapped)
total_unmapped <- sum(SSU$unmapped)

mapped_SSU <- sum(SSU$mapped)
mapped_LSU <- sum(LSU$mapped)


# match 
match = fread('/mnt/data3/yiyonghao/MicroRNA/code_upload/4.rRNAcalculated/summary.tsv')
unmapped_matched   <- match$matched_rrna_qseqids
unmapped_unmatched <- match$unmatched_qseqids


# blast
library(tidyr)
blast = fread('/mnt/data3/yiyonghao/MicroRNA/code_upload/4.rRNAcalculated/details.tsv',header = TRUE)
max_len = max(sapply(strsplit(blast$matched_terms,';'),length))
blast = blast %>% 
    separate(matched_terms,into = paste0('col',1:max_len),sep = ';',fill = 'right') %>% 
    select(-file)

blast_long = blast %>% 
    pivot_longer(cols = -qseqid,names_to = 'index',values_to = 'value'
    ) %>% na.omit()
table(blast_long$value)
tbl <- c(`16S`=9349, `18S`=1589, `23S`=15521, `28S`=1810, `5.8S`=1562, `5S`=12865)


# draw 
rrna_nested_donut_nolabel <- function(
  total_mapped,
  total_unmapped,
  mapped_SSU,
  mapped_LSU,
  unmapped_matched,
  unmapped_unmatched,
  unmapped_matched_breakdown,   # 命名向量，如 c("16S"=..., "18S"=..., ...)
  inner_r0 = 0.70,
  inner_r  = 1.20,
  outer_gap = 0.05,
  outer_thickness = 0.65,
  palette = c(
    mapped   = "#4C78A8",
    unmapped = "#E45756",
    SSU      = "#72B7B2",
    LSU      = "#F58518",
    unmapped_unmatched = "#A0A0A0"
  ),
  title = "rRNA overview",
  save_path = NULL, width = 7, height = 7,
  return_data = FALSE
){
  total_all <- total_mapped + total_unmapped
  if (total_all <= 0) stop("总量为 0，无法绘图。")

  o_r0 <- inner_r + outer_gap
  o_r  <- o_r0 + outer_thickness

  # 内圈
  frac_mapped <- total_mapped / total_all
  inner_df <- tibble::tibble(
    level=1L, parent=NA_character_,
    label=c("mapped","unmapped"),
    value=c(total_mapped,total_unmapped),
    start=c(0, 2*pi*frac_mapped),
    end=c(2*pi*frac_mapped, 2*pi),
    r0=inner_r0, r=inner_r
  )

  m_start <- inner_df$start[inner_df$label=="mapped"];   m_end <- inner_df$end[inner_df$label=="mapped"]
  u_start <- inner_df$start[inner_df$label=="unmapped"]; u_end <- inner_df$end[inner_df$label=="unmapped"]
  m_len <- m_end - m_start; u_len <- u_end - u_start

  # mapped 外圈
  mapped_children <- tibble::tibble(
    level=2L,parent="mapped",
    label=c("SSU","LSU"),
    value=c(mapped_SSU,mapped_LSU)
  ) %>%
    mutate(frac=value/sum(value),
           end=m_start+m_len*cumsum(frac),
           start=lag(end,default=m_start),
           r0=o_r0,r=o_r)

  # unmapped 外圈
  unmapped_split <- tibble::tibble(
    level=2L,parent="unmapped",
    label=c("unmapped_unmatched","unmapped_matched_total"),
    value=c(unmapped_unmatched,unmapped_matched)
  ) %>%
    mutate(frac=value/sum(value),
           end=u_start+u_len*cumsum(frac),
           start=lag(end,default=u_start),
           r0=o_r0,r=o_r)

  # matched 内细分
  detail_df <- tibble()
  if (!is.null(unmapped_matched_breakdown) && sum(unmapped_matched_breakdown)>0) {
    rel <- unmapped_matched_breakdown/sum(unmapped_matched_breakdown)
    scaled <- rel*unmapped_matched
    ms <- unmapped_split %>% filter(label=="unmapped_matched_total")
    ms_len <- ms$end-ms$start
    detail_df <- tibble::tibble(
      level=2L,parent="unmapped",
      label=names(scaled), value=as.numeric(scaled)
    ) %>%
      mutate(frac=value/sum(value),
             end=ms$start+ms_len*cumsum(frac),
             start=lag(end,default=ms$start),
             r0=o_r0,r=o_r)
  }

  arcs <- bind_rows(inner_df, mapped_children,
                    unmapped_split %>% filter(label=="unmapped_unmatched"),
                    detail_df)

  # 调色
  cats <- unique(arcs$label)
  pal <- palette
  miss <- setdiff(cats, names(pal))
  if (length(miss)>0) pal <- c(pal, setNames(hcl.colors(length(miss),"Dark 3"), miss))

  # 作图：仅环图+图例，无文字标签
  p <- ggplot(arcs) +
    geom_arc_bar(aes(x0=0,y0=0,r0=r0,r=r,start=start,end=end,fill=label),
                 color="white",size=0.4) +
    coord_fixed(clip="off") +
    theme_void() +
    scale_fill_manual(values=pal) +
    guides(fill=guide_legend(title=NULL)) +
    ggtitle(title)

  if (!is.null(save_path)) {
    ggsave(save_path, plot=p, width=width, height=height, dpi=300)
  }
  if (isTRUE(return_data)) return(list(plot=p,data=arcs))
  return(p)
}
p <- rrna_nested_donut_nolabel(
  total_mapped, total_unmapped,
  mapped_SSU, mapped_LSU,
  unmapped_matched, unmapped_unmatched,
  unmapped_matched_breakdown = tbl,
  title = ""
)
p

ggsave(plot = p,
        filename = '/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/rRNA_ratio.pdf')
mapped_ratio = total_mapped / (total_mapped + total_unmapped)
mapped_ratio
unmapped_ratio = total_unmapped / (total_mapped + total_unmapped)
unmapped_ratio

LSU_ratio = mapped_LSU / total_mapped
LSU_ratio
SSU_ratio = mapped_SSU / total_mapped
SSU_ratio

match_sum = 9439+1589+15521+1810+1562+12865
for(i in 1:length(tbl)){
  ratio = ((tbl[i] / match_sum) * unmapped_matched) / (unmapped_matched + unmapped_unmatched)
  print(paste0(names(tbl[i]),'_',ratio))
}

unmatched = unmapped_unmatched / (unmapped_unmatched+unmapped_matched)
unmatched
