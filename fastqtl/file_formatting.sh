2. Troubleshooting fastQTL - ERROR: Genotype data does not overlap with phenotype data, check your files!

	***** identified the the names in the bed and sample files did not match the VCF *****

		***1] Replace the names in the AA_eqtl.txt file to match the VCF 

gzip -cd Whole_Blood_Analysis.v6p.normalized.expression.bed.gz | head -1 | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX'

gzip -cd Whole_Blood_Analysis.v6p.normalized.expression.bed.gz | head -1 | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX' | wc -l

gzip -cd Whole_Blood_Analysis.v6p.normalized.expression.bed.gz | head -1 | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX' > bed_samples.txt

gzip -cd Whole_Blood_Analysis.v6p.normalized.expression.bed.gz | head -1 | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX' > bed_samples.txt

gzip -cd GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz | head -10000 | grep -v '##' | grep '#' | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX' | head

grep 'GTEX-111YS' vcf_samples.txt 

gzip -cd GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz | head -10000 | grep -v '##' | grep '#' | tr '\t' '\n' | grep -vE '^$' | grep 'GTEX' > vcf_samples.txt

wc -l vcf_samples.txt

head vcf_samples.txt

head bed_samples.txt 

grep 'GTEX-111YS' vcf_samples.txt 

yocelyn@instance-1:~/fastQTL$ head Gtex_AA.txt
GTEX-11DXZ
GTEX-11DYG
GTEX-11NV4
GTEX-11P7K
GTEX-11PRG
GTEX-11TT1
GTEX-11WQC
GTEX-12126
GTEX-12WSD
GTEX-131XE

yocelyn@instance-1:~/fastQTL$ grep -f Gtex_AA.txt vcf_samples.txt | wc -l
50

yocelyn@instance-1:~/fastQTL$ grep -f Gtex_AA.txt vcf_samples.txt 
GTEX-SNOS-0004-SM-325AE
GTEX-VJYA-0003-SM-3G3S9
GTEX-R55D-0926-SM-325AK
GTEX-WCDI-0001-SM-3PZ7A
GTEX-QLQW-0001-SM-2MRID
GTEX-U412-0001-SM-3G3S1
GTEX-WXYG-0001-SM-3PZ82
GTEX-S341-0001-SM-2VCTY
GTEX-Q2AH-0004-SM-2IJGI
GTEX-TML8-0004-SM-325A2
GTEX-UJMC-0002-SM-3G3ST
GTEX-POMQ-0003-SM-2H2X7
GTEX-RU1J-0002-SM-2RXF5
GTEX-T2YK-0003-SM-325AB
GTEX-W5X1-0002-SM-3G3S7
GTEX-TKQ1-0002-SM-3G3SW
GTEX-SSA3-0001-SM-325A5
GTEX-X585-0003-SM-3PZ7U
GTEX-R53T-0003-SM-2RXEX
GTEX-OHPK-0003-SM-2H2WX
GTEX-X3Y1-0002-SM-3PZ85
GTEX-VUSH-0002-SM-3PZ7N
GTEX-P4QS-0003-SM-2H2X1
GTEX-UJHI-0001-SM-3L287
GTEX-NFK9-0004-SM-2BWYG
GTEX-Y114-0004-SM-4YUWP
GTEX-11TT1-0001-SM-5DWT7
GTEX-13OVI-0003-SM-5HJQK
GTEX-XYKS-0003-SM-4YUVD
GTEX-ZQUD-0002-SM-4YUWT
GTEX-11DXZ-0003-SM-5DWSR
GTEX-YFCO-0001-SM-4YUWM
GTEX-13FXS-0004-SM-5HJPR
GTEX-13VXT-0003-SM-5HJR7
GTEX-13JUV-0001-SM-5HJPQ
GTEX-12WSD-0002-SM-5DWRK
GTEX-ZVP2-0002-SM-4YUXK
GTEX-ZV7C-0002-SM-4YUX7
GTEX-11NV4-0001-SM-5DWSY
GTEX-XMD3-0003-SM-4YUX3
GTEX-11WQC-0003-SM-5DWTG
GTEX-11P7K-0002-SM-5DWTW
GTEX-11PRG-0004-SM-5DWRB
GTEX-11DYG-0004-SM-5DWSF
GTEX-ZPU1-0002-SM-4YUW5
GTEX-12126-0001-SM-5DWTE
GTEX-ZTX8-0004-SM-4YUXI
GTEX-13O61-0003-SM-5HJRK
GTEX-Y3IK-0004-SM-4YUX2
GTEX-ZVTK-0002-SM-4YUW7


yocelyn@instance-1:~/fastQTL$ grep -f Gtex_EUR.txt vcf_samples.txt 
GTEX-P4PP-0004-SM-2H2WW
GTEX-PW2O-0004-SM-2IJGJ
GTEX-U3ZH-0003-SM-3G3RZ
GTEX-PLZ4-0003-SM-2IJGB
GTEX-WOFL-0003-SM-3PZ7H
GTEX-QMRM-0002-SM-2MRIF
GTEX-Q2AG-0003-SM-2IJGG
GTEX-QDT8-0926-SM-325A6
GTEX-QV31-0003-SM-2MRI2
GTEX-PLZ6-0003-SM-2IJGD
GTEX-TMMY-0003-SM-3G3SV
GTEX-WYJK-0001-SM-3PZ7Z
GTEX-Q2AI-0003-SM-2IJGK
GTEX-XBED-0004-SM-3PZ7G
GTEX-RUSQ-0001-SM-2RXER
GTEX-U8XE-0004-SM-3G3SU

yocelyn@instance-1:~/fastQTL$ grep -f Gtex_EUR.txt vcf_samples.txt | wc -l
242
yocelyn@instance-1:~/fastQTL$ grep -f Gtex_EUR.txt vcf_samples.txt > Gtex_EUR_long.txt

	***2] Replace the AA and Eur names in the BED files with format that matches VCF  IN Replace

bed = read.table('Whole_Blood_Analysis.v6p.normalized.expression.bed', header = TRUE, comment.char = '', as.is = TRUE)

bed[1:10, 1:5]
   X.chr  start    end           gene_id GTEX.111YS
1      1  11868  11869 ENSG00000223972.4 -0.2655829
2      1  29552  29553 ENSG00000227232.4 -1.1850740
3      1 129222 129223 ENSG00000238009.2  0.1671441
4      1 131024 131025 ENSG00000233750.3  0.2809339
5      1 135894 135895 ENSG00000268903.1  0.4966030
6      1 139378 139379 ENSG00000237683.5  0.5134028


> colnames(bed) <- c('#chr', colnames(bed)[2:ncol(bed)])
> bed[1:10, 1:5]
   #chr  start    end           gene_id GTEX.111YS
1     1  11868  11869 ENSG00000223972.4 -0.2655829
2     1  29552  29553 ENSG00000227232.4 -1.1850740
3     1 129222 129223 ENSG00000238009.2  0.1671441


> vcfnames = readLines('vcf_samples.txt')
> vcfnames
  [1] "GTEX-P4PP-0004-SM-2H2WW"    "GTEX-PWO3-0004-SM-2IJGP"   
  [3] "GTEX-SNOS-0004-SM-325AE"    "GTEX-PW2O-0004-SM-2IJGJ"   
  [5] "GTEX-RU72-0004-SM-2RXF6"    "GTEX-VJYA-0003-SM-3G3S9"   
  [7] "GTEX-U3ZH-0003-SM-3G3RZ"    "GTEX-PLZ4-0003-SM-2IJGB"   
  [9] "GTEX-P44G-0926-SM-325AM"    "GTEX-WOFL-0003-SM-3PZ7H"   
 [11] "GTEX-R55D-0926-SM-325AK"    "GTEX-QMRM-0002-SM-2MRIF"   
 [13] "GTEX-Q2AG-0003-SM-2IJGG"    "GTEX-WCDI-0001-SM-3PZ7A"  


 > length(vcfnames)
[1] 450

pmatch(
x=              table=          nomatch=        duplicates.ok=  
> pmatch(colnames(bed), vcfnames)
  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA


  > head(colnames(bed)
+ )
[1] "#chr"       "start"      "end"        "gene_id"    "GTEX.111YS"
[6] "GTEX.1122O"
> head(gsub('.', '-', colnames(bed), fixed = TRUE))


> na.omit(pmatch(colnames(bed), vcfnames))
  [1] 277 242 240 234 379 430 398 410 402 224 260 239 406 397 252 387 437 328
 [19] 241 216 220 380 189 195 377 447 401 245 282 403 191 442 307 207 420 390
 [37] 350 320 206 308 291 185 427 445 232 399 337 374 318 342 444 435 198 225
 [55] 192 326 349 376 419 365 312 249 208 288 276 186 268 247 439 370 188 273
 [73] 296 238 385 363 274 259 270 281 343 318 309 339 434 251 389 201 369 432
 [91] 433 405 280 327 375 264 407 261 336 341 226 109 180 185 267  58  75  85
[109] 183 163 142 188  91 153  56  19 167 103  41  57 152 131  50  98 181  99
[127]   1 150 119 102   8 144  18  87 139 137  71   4 181 129   2 114 157  13
[145]  64  21  42 151  15  70  26  51 141  29  34  12 182  47 187  74 189 138
[163]  96  11  61 149  36 161 115 112  89   5 182 110 148 169 101  59 156 160
[181]  30  39  81 186  54 162 136  32 135  31  73   3 122  95  90 120  62  72
[199] 130 158 121 159  68  18  69 125 196   7  80  92  28  24 184  86 100  44
[218] 164  83  37  93   6  27 166  46 118  14  60 146  49  38  53 124  43 123
[235]  94 147 107  33  10 126  82  65 132 105 165  55  76  20 127 113 116 104
[253] 143  66 106  35 111 134 128  25 140 324 394  45 118 450  22 299 193 348
[271] 329 382 351 289 313 354 227 256 353 381 292 253 310 204 315 183 378 443
[289] 266 194 316 431 357 388 393 411 362 325 231 446 422 286 441 344 255 356
[307] 218 416 287 212 205 230 384 355 400 368 209 418 332 213 243 184 283 391
[325] 429 360 414 366 436 364 367 295 190 449 237 383 392 263

attr(,"na.action")
[1] 1 2 3 4
attr(,"class")
[1] "omit"
> head(vcfnames[na.omit(pmatch(colnames(bed), vcfnames))])
[1] "GTEX-111YS-0004-SM-5DWSS" "GTEX-1122O-0004-SM-5DWSJ"
[3] "GTEX-1128S-0001-SM-5DWTV" "GTEX-113IC-0001-SM-5DWSI"
[5] "GTEX-113JC-0001-SM-5DWSH" "GTEX-118XS-0001-SM-5DWSE"
> head(colnames(bed)[4:ncol(bed)])
[1] "gene_id"    "GTEX-111YS" "GTEX-1122O" "GTEX-1128S" "GTEX-113IC"
[6] "GTEX-113JC"
> head(colnames(bed)[5:ncol(bed)])
[1] "GTEX-111YS" "GTEX-1122O" "GTEX-1128S" "GTEX-113IC" "GTEX-113JC"
[6] "GTEX-118XS"
> long.names = vcfnames[na.omit(pmatch(colnames(bed), vcfnames))]
> length(long.names)
[1] 338
> ncol(bed) - 4
[1] 338
> colnames(bed)[1:4]
[1] "#chr"    "start"   "end"     "gene_id"
> c(colnames(bed)[1:4], long.names)
  [1] "#chr"                       "start"                     
  [3] "end"                        "gene_id"                   
  [5] "GTEX-111YS-0004-SM-5DWSS"   "GTEX-1122O-0004-SM-5DWSJ"  
  [7] "GTEX-1128S-0001-SM-5DWTV"   "GTEX-113IC-0001-SM-5DWSI"

> colnames(bed) = c(colnames(bed)[1:4], long.names)
> bed[1:10, 1:5]
> write.table(x = bed, file = 'Whole_Blood_Analysis.v6p.normalized.expression.long_names.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

bgzip Whole_Blood_Analysis.v6p.normalized.expression.long_names.bed
tabix -p bed Whole_Blood_Analysis.v6p.normalized.expression.long_names.bed.gz 