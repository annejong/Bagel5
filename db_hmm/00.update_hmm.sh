# PFAM: wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# InterPro to HMM: https://www.ebi.ac.uk/interpro/entry/InterPro/#table   ==> export all
# TIGR: wget       https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv


# db_hmm folder
db_hmm_dir=/data/bagel5/db_hmm


# Inherited from bagel4
cat $db_hmm_dir/BAGEL4/*.hmm > $db_hmm_dir/bagel4_bacteriocins.hmm
hmmpress -f $db_hmm_dir/bagel4_bacteriocins.hmm

# LinL for AOI discovery
cp $db_hmm_dir/BAGEL4/LinL.hmm $db_hmm_dir/LinL.hmm 
hmmpress -f $db_hmm_dir/LinL.hmm



# Get PFAM models
function Get_PFAM {
python3 /data/bagel5/scripts/extract_hmm_from_pfamdb.py -pfamdb /data/databases/PFAM/Pfam-A.hmm -outdir $db_hmm_dir/PFAM -acc $1
}

Get_PFAM PF00005
Get_PFAM PF00067
Get_PFAM PF00069
Get_PFAM PF00069
Get_PFAM PF00072
Get_PFAM PF00082
Get_PFAM PF00107
Get_PFAM PF00296
Get_PFAM PF00310
Get_PFAM PF00486
Get_PFAM PF00512
Get_PFAM PF00535
Get_PFAM PF00535
Get_PFAM PF00733
Get_PFAM PF00733
Get_PFAM PF01320
Get_PFAM PF01494
Get_PFAM PF01636
Get_PFAM PF01636
Get_PFAM PF01721
Get_PFAM PF02441
Get_PFAM PF02624
Get_PFAM PF02624
Get_PFAM PF02794
Get_PFAM PF03047
Get_PFAM PF03412
Get_PFAM PF03412
Get_PFAM PF03857
Get_PFAM PF04055
Get_PFAM PF04055
Get_PFAM PF04369
Get_PFAM PF04738
Get_PFAM PF04738
Get_PFAM PF05147
Get_PFAM PF05147
Get_PFAM PF05402
Get_PFAM PF07475
Get_PFAM PF07812
Get_PFAM PF07812
Get_PFAM PF08109
Get_PFAM PF08129
Get_PFAM PF08951
Get_PFAM PF09221
Get_PFAM PF09683
Get_PFAM PF10439
Get_PFAM PF11420
Get_PFAM PF11632
Get_PFAM PF11758
Get_PFAM PF12173
Get_PFAM PF12679
Get_PFAM PF12730
Get_PFAM PF13186
Get_PFAM PF13186
Get_PFAM PF13437
Get_PFAM PF13471
Get_PFAM PF13471
Get_PFAM PF13537
Get_PFAM PF13537
Get_PFAM PF13581
Get_PFAM PF13632
Get_PFAM PF14028
Get_PFAM PF14028 
Get_PFAM PF16942
Get_PFAM PF17312
Get_PFAM PF17508
Get_PFAM PF17914
Get_PFAM PF17914
Get_PFAM PF19409

# Get TIGR models
function Get_TIGR {
wget https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/$1.HMM  -P $db_hmm_dir/TIGR
}

Get_TIGR TIGR00702.1
Get_TIGR TIGR03605.1
Get_TIGR TIGR03895.1
Get_TIGR TIGR03897.1
Get_TIGR TIGR03975.1
Get_TIGR TIGR04078.1
Get_TIGR TIGR04184.1
Get_TIGR TIGR04185.1
Get_TIGR TIGR04445.1
Get_TIGR TIGR04517.1
Get_TIGR TIGR03732.1	
Get_TIGR TIGR03733.1	
Get_TIGR TIGR03740.1	
Get_TIGR TIGR04085.1	
Get_TIGR TIGR03891.1	
Get_TIGR TIGR03892.1
Get_TIGR TIGR01847.1
Get_TIGR TIGR03898.1
Get_TIGR TIGR03973.1
Get_TIGR TIGR04079.1
Get_TIGR TIGR04149.1
Get_TIGR TIGR04363.1
Get_TIGR TIGR03731.1
Get_TIGR TIGR03893.1
Get_TIGR NF037999.1
Get_TIGR NF040663.1
Get_TIGR NF000539.0
Get_TIGR NF033457.1
Get_TIGR TIGR01653.1
Get_TIGR TIGR03601.1
Get_TIGR TIGR03795.1
Get_TIGR TIGR03892.1
Get_TIGR TIGR04067.1
Get_TIGR NF013850.1
Get_TIGR NF015035.4
Get_TIGR NF016272.4
Get_TIGR NF019720.4
Get_TIGR NF019739.4
Get_TIGR NF020782.4
Get_TIGR NF021219.4
Get_TIGR NF021923.4
Get_TIGR NF022861.4
Get_TIGR NF023063.4
Get_TIGR NF023186.4
Get_TIGR NF023595.4
Get_TIGR NF028252.4
Get_TIGR NF033222.1
Get_TIGR NF033384.1
Get_TIGR NF033385.1
Get_TIGR NF033400.3
Get_TIGR NF033557.1
Get_TIGR NF033806.2
Get_TIGR NF033837.1
Get_TIGR NF033881.1
Get_TIGR NF035954.1
Get_TIGR NF036357.4
Get_TIGR NF036612.4
Get_TIGR NF038034.1
Get_TIGR NF038035.1
Get_TIGR NF038163.1
Get_TIGR NF038164.1
Get_TIGR NF038165.1
Get_TIGR NF038168.1
Get_TIGR NF039937.3
Get_TIGR NF033433.1	


############################################# 
#    HMMs for AOI discovery                 #
#############################################


# PFAMs for AOI discovery
cat $db_hmm_dir/PFAM/PF00069.hmm  >$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00310.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00535.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00733.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF01636.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF02624.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF03412.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF04055.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF05147.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF07812.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13186.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13471.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13537.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13632.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF17914.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF04738.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
cat $db_hmm_dir/PFAM/PF14028.hmm >>$db_hmm_dir/AOI_discovery_PFAM.hmm
hmmpress -f $db_hmm_dir/AOI_discovery_PFAM.hmm

# TIGRs for AOI discovery
cat $db_hmm_dir/TIGR/TIGR00702.1.HMM   >$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03605.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03895.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03897.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03975.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04078.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04184.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04185.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04445.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04517.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03732.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03733.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03740.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04085.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03891.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03892.1.HMM  >>$db_hmm_dir/AOI_discovery_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033433.1.HMM   >>$db_hmm_dir/AOI_discovery_TIGR.hmm
hmmpress -f $db_hmm_dir/AOI_discovery_TIGR.hmm




############################################# 
#    HMMs for AOI annotation                #
#############################################



# PFAMs for AOI Annotation Core Peptides

cat $db_hmm_dir/PFAM/PF01721.hmm  >$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF03047.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF04369.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF08109.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF08129.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF09221.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF09683.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF10439.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF11420.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF11632.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF11758.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF12173.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF16942.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF17312.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF17508.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
cat $db_hmm_dir/PFAM/PF19409.hmm >>$db_hmm_dir/AOI_Annotation_Core_PFAM.hmm
hmmpress -f $db_hmm_dir/AOI_Annotation_Core_PFAM.hmm


# PFAMs for AOI Annotation Context
cat $db_hmm_dir/PFAM/PF00069.hmm  >$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00310.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00535.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00733.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF01636.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF02624.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF03412.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF04055.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF05147.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF07812.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13186.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13471.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13537.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13632.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF17914.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF04738.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF14028.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00005.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00067.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00072.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00082.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00107.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00296.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00486.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF00512.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF01320.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF01494.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF02441.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF02794.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF03857.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF05402.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF07475.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF08951.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF12679.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF12730.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13437.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
cat $db_hmm_dir/PFAM/PF13581.hmm >>$db_hmm_dir/AOI_Annotation_Context_PFAM.hmm
hmmpress -f $db_hmm_dir/AOI_Annotation_Context_PFAM.hmm

# TIGRs for AOI Annotation Context
cat $db_hmm_dir/TIGR/TIGR00702.1.HMM   >$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03605.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03895.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03897.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03975.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04078.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04184.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04185.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04445.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04517.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03732.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03733.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03740.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04085.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03891.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033433.1.HMM   >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR00702.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04184.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04185.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04445.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04517.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03732.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03733.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03740.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03891.1.HMM  >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033433.1.HMM   >>$db_hmm_dir/AOI_Annotation_Context_TIGR.hmm
hmmpress -f $db_hmm_dir/AOI_Annotation_Context_TIGR.hmm

# TIGRs for AOI Annotation Core Peptide
cat $db_hmm_dir/TIGR/NF000539.0.HMM    >$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033222.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033384.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033400.3.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033457.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033557.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033806.2.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033837.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF033881.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF035954.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF037999.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF038034.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF038035.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF038163.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF038164.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF038168.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/NF040663.1.HMM   >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR01653.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR01847.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03731.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03795.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03892.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03892.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03893.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03898.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR03973.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04067.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04079.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
cat $db_hmm_dir/TIGR/TIGR04149.1.HMM  >>$db_hmm_dir/AOI_Annotation_Core_TIGR.hmm
hmmpress -f $db_hmm_dir/AOI_Annotation_Core_TIGR.hmm


# BAGEL4 motifs for AOI Annotation Core Peptide
cat $db_hmm_dir/BAGEL4/bacteriocinII.hmm >$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/ggmotif.hmm       >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/leaderLanBC.hmm   >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/leaderLanM.hmm    >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/LE-DUF.hmm        >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/LE-LAC481.hmm     >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/LE-LanBC.hmm      >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/LE-MER+2PEP.hmm   >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/LinL.hmm          >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-2PEPA.hmm      >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-2PEPB.hmm      >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-DUF.hmm        >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-EPI.hmm        >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-LAC481.hmm     >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-NIS+EPI.hmm    >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
cat $db_hmm_dir/BAGEL4/MA-NIS.hmm        >>$db_hmm_dir/AOI_Annotation_Core_BAGEL4.hmm
