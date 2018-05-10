1. Extract Select Samples 
expr = read.table('Whole_Blood_Analysis.v6p.normalized.expression.bed', header = TRUE, sep = '\t')
list.files('~/Downloads', pattern = 'Gtex')
afr = readLines('~/Downloads/Gtex_AA.txt')
match(afr, colnames(expr))
colnames(expr)
expr = read.table('Whole_Blood_Analysis.v6p.normalized.expression.bed', header = TRUE, sep = '\t', comm)
expr = read.table('Whole_Blood_Analysis.v6p.normalized.expression.bed', header = TRUE, sep = '\t', comment.char = '')
> colnames(expr)
afr = gsub('-', '.', afr)
match(afr, colnames(expr))
c(1:4, na.omit(match(afr, colnames(expr))))
> afr.use = c(1:4, na.omit(match(afr, colnames(expr))))
> colnames(expr)[1:4]
> afr_expr = expr[, afr.use]
> dim(afr_expr)
[1] 23152    52
> length(na.omit(match(afr, colnames(expr))))
[1] 48
> head(afr_expr)
> gsub('.', '-', colnames(afr_expr), fixed = TRUE)
> colnames(afr_expr) = gsub('.', '-', colnames(afr_expr), fixed = TRUE)
> colnames(afr_expr)
> colnames(afr_expr)[2:ncol(afr_expr)]
> c('#chr', colnames(afr_expr)[2:ncol(afr_expr)])
> colnames(afr_expr) = c('#chr', colnames(afr_expr)[2:ncol(afr_expr)])
> colnames(afr_expr)[2:ncol(afr_expr)]
> colnames(afr_expr)
 [1] "#chr"       "start"      "end"        "gene_id"    "GTEX-11DXZ"
 [6] "GTEX-11DYG" "GTEX-11NV4" "GTEX-11P7K" "GTEX-11PRG" "GTEX-11TT1"
[11] "GTEX-11WQC" "GTEX-12126" "GTEX-12WSD" "GTEX-13FXS" "GTEX-13JUV"
[16] "GTEX-13O61" "GTEX-13OVI" "GTEX-13VXT" "GTEX-NFK9"  "GTEX-OHPK" 
[21] "GTEX-POMQ"  "GTEX-Q2AH"  "GTEX-QLQW"  "GTEX-R53T"  "GTEX-R55D" 
[26] "GTEX-RU1J"  "GTEX-S341"  "GTEX-SNOS"  "GTEX-SSA3"  "GTEX-T2YK" 
[31] "GTEX-TKQ1"  "GTEX-TML8"  "GTEX-UJHI"  "GTEX-UJMC"  "GTEX-VJYA" 
[36] "GTEX-VUSH"  "GTEX-W5X1"  "GTEX-WCDI"  "GTEX-WXYG"  "GTEX-X3Y1" 
[41] "GTEX-X585"  "GTEX-XMD3"  "GTEX-XYKS"  "GTEX-Y114"  "GTEX-Y3IK" 
[46] "GTEX-YFCO"  "GTEX-ZPU1"  "GTEX-ZQUD"  "GTEX-ZTX8"  "GTEX-ZV7C" 
[51] "GTEX-ZVP2"  "GTEX-ZVTK" 