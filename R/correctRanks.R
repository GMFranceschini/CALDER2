#' correctCalderRank
#'
#' @description
#' Implements a correction of continuous ranks.
#' Writes a new .bed file in the same location as the original one with the fixed rank.
#' The new file has the exact same structure and column names of the original, maximizing compatibility with
#' the rest of existent code base.
#' 
#' @param calder_call the path to the CALDER call output (bed file)
#'
#' @return Nothing
#' @export
correctCalderRank = function(calder_call){
  
  df = fread(
    calder_call,
    sep = "\t",
    quote = "",
    header = F,
    data.table = F,
    stringsAsFactors = FALSE
  )
  
  if (ncol(df) == 9){
    
    coln = c(
      "chr",
      "start",
      "end",
      "comp",
      "comp_rank",
      "empty",
      "start2",
      "end2",
      "color",
      "boh",
      "bin_comp"
    )} else {
      coln = c(
        "chr",
        "start",
        "end",
        "comp",
        "comp_rank",
        "empty",
        "start2",
        "end2",
        "color",
        "boh",
        "bin_comp"
      )
    }
  
  colnames(df) = coln
  
  df = df %>% mutate(l8 = substr(comp, 1,5))
  
  mm = df %>% 
    group_by(l8, chr) %>% 
    summarise(old_min = min(comp_rank), old_max = max(comp_rank))
  
  rr = data.frame(
    l8 = c("A.1.1", "A.1.2", "A.2.1", "A.2.2", "B.1.1", "B.1.2", "B.2.1", "B.2.2"),
    new_min = rev(seq(0, 7/8, 1/8)),
    new_max = rev(seq(1/8, 1, 1/8)))
  
  df = df %>% inner_join(rr, by = "l8") %>% 
    inner_join(mm, by = c("chr", "l8")) %>% 
    mutate(new_CompRank = new_min + ((comp_rank-old_min)/(old_max-old_min))*(new_max-new_min)) 
  
  df_new = df %>% 
    mutate(comp_rank = new_CompRank) %>% 
    select(all_of(coln)) %>% 
    arrange(chr, start)
  
  name_new = gsub(".bed", "_corrRank.bed", calder_call)
  
  message("Writing corrected file...")
  
  fwrite(
    df_new,
    file = name_new,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  
}
