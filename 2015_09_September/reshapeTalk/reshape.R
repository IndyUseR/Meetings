### R code from vignette source 'reshape.Rnw'

###################################################
### code chunk number 1: load_data
###################################################
library(reshape)
library(Hmisc)

load("mtx.RData")
################################################################################
# [1] "patid"               "visitid"             "dob"                
# [4] "trt_age"             "sex"                 "height"             
# [7] "weight"              "ibw"                 "adjust_bw"          
#[10] "bsa"                 "dosing_bsa"          "dx"                 
#[13] "pleural_effusion"    "ascites"             "mtxlvl_0"           
#[16] "mtx_mgm2"            "infuse_timeh"        "infuse_timecode"    
#[19] "mtx_startdt"         "mtx_starttm"         "leu_init_dose"      
#[22] "leu_init_freq"       "leu_route"           "keppra"             
#[25] "keppra_dose"         "keppra_freq"         "num_inter_meds"     
#[28] "hydrate_fluids"      "hydrate_fluids_rate" "mtx_startdt_base"   
#[31] "mtx_starttm_base"    "urine_0"             "scr_0"              
#[34] "crclwt_0"            "crcl_0"              "ast_0"              
#[37] "alt_0"               "labdt_24"            "labtm_24"           
#[40] "urine_24"            "scr_24"              "crcl_24"            
#[43] "mtxlvl_24"           "delaymtx_24"         "labdt_48"           
#[46] "labtm_48"            "urine_48"            "scr_48"             
#[49] "crcl_48"             "mtxlvl_48"           "delaymtx_48"        
#[52] "labdt_72"            "labtm_72"            "urine_72"           
#[55] "scr_72"              "crcl_72"             "mtxlvl_72"          
#[58] "delaymtx_72"         "labdtadd_1"          "labtmadd_1"         
#[61] "urineadd_1"          "scradd_1"            "crcladd_1"          
#[64] "mtxlvladd_1"         "labdtadd_2"          "labtmadd_2"         
#[67] "urineadd_2"          "scradd_2"            "crcladd_2"          
#[70] "mtxlvladd_2"         "labdtadd_3"          "labtmadd_3"         
#[73] "urineadd_3"          "scradd_3"            "crcladd_3"          
#[76] "mtxlvladd_3"         "delaymtx_any"        "nontoxic_mtx_tm"    
#[79] "urine_add3"          "mtx_lvl_add1"        "mtx_lvl_24h"        
#[82] "delaymtx2_any"      
################################################################################
mtx_both <- mtx_dt[,c("patid","visitid","keppra","trt_age","mtxlvl_24","scr_24","mtxlvl_48","scr_48","mtxlvl_72","scr_72")]
names(mtx_both)[5:10] <- c("mtxlvl.24","scr.24","mtxlvl.48","scr.48","mtxlvl.72","scr.72")
mtx.mltf <- melt(mtx_both, id=c("patid","visitid","keppra","trt_age"),  variable_name="lab.var")
mtx.mltf.srt <- with(mtx.mltf, mtx.mltf[order(patid,visitid),])
mtx.mltf.lab <- cbind(mtx.mltf.srt, colsplit(mtx.mltf.srt$lab.var, names=c("lab","time"), split="\\."))

mtx.mtxlvl <- mtx_dt[,c("patid","visitid","keppra","trt_age","mtxlvl_24","mtxlvl_48","mtxlvl_72")]
names(mtx.mtxlvl)[5:7] <- c("mtxlvl.24","mtxlvl.48","mtxlvl.72")
mtx.mlt1 <- melt(mtx.mtxlvl, id=c("patid","visitid","keppra","trt_age"), measured=c("mtxlvl.24","mtxlvl.48","mtxlvl.72"), variable_name="lab.var")
mtx.mlt1.srt <- with(mtx.mlt1, mtx.mlt1[order(patid,visitid),])
mtx.mlt1.lab <- cbind(mtx.mlt1.srt, colsplit(mtx.mlt1.srt$lab.var, names=c("lab","time"), split="\\."))


###################################################
### code chunk number 2: mtxraw
###################################################
latex(mtx.mtxlvl[1:5,c(1:3,5:7)],file="",rowname=NULL)


###################################################
### code chunk number 3: mtxmelt
###################################################
latex(mtx.mlt1.srt[1:6,c(1:3,5:6)], file="", rowname=NULL)


###################################################
### code chunk number 4: mlttrt
###################################################
latex(mtx.mlt1.lab[1:3,c(1:3,6:8)], file="", rowname=NULL)


###################################################
### code chunk number 5: meas2
###################################################
latex(head(mtx.mltf.lab[1:6,c(1:3,5:6)]), file="", rowname=NULL)


###################################################
### code chunk number 6: mltmeas2
###################################################
latex(head(mtx.mltf.lab[1:6,c(1:3,5:8)]), file="", rowname=NULL)


###################################################
### code chunk number 7: castit
###################################################
latex(head(cast(mtx.mlt1.srt, patid+visitid~lab.var)), file="", rowname=NULL)


###################################################
### code chunk number 8: castit2
###################################################
mtx.castf.lab <- cast(mtx.mltf.lab, patid+visitid+keppra+time~lab)
latex(head(mtx.castf.lab), file="", rowname=NULL)


###################################################
### code chunk number 9: castag
###################################################
latex(cast(mtx.mltf.lab, keppra~lab, mean, na.rm=TRUE), file="", rowname=NULL, digits=2)
#cast(mtx.mltf.lab, keppra+time~lab, mean, na.rm=TRUE, margins="keppra")
#cast(mtx.mltf.lab, keppra~lab, mean, na.rm=TRUE, margins="grand_col")
#cast(mtx.mltf.lab, keppra~lab, mean, na.rm=TRUE, margins="grand_row")
#cast(mtx.mltf.lab, keppra~lab, mean, na.rm=TRUE, margins=TRUE)


###################################################
### code chunk number 10: castag2
###################################################
latex(cast(mtx.mltf.lab, lab~., quantile, c(0.25,0.5,0.75),na.rm=TRUE), file="", rowname=NULL)


###################################################
### code chunk number 11: castag3
###################################################
latex(cast(mtx.mltf.lab, time+keppra~lab, mean, na.rm=TRUE), file="", rowname=NULL, digits=2)


###################################################
### code chunk number 12: castag4
###################################################
latex(cast(mtx.mltf.lab, time+keppra~lab, mean, na.rm=TRUE, margins="time"), file="", rowname=NULL, digits=2)


