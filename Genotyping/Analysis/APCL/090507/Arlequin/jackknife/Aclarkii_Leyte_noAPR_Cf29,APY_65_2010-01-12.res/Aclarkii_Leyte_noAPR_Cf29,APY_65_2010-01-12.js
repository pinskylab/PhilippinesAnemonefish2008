
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Documents%20and%20Settings/Administrator/My%20Documents/PopGen%20Programs/Arlequin311/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 12/01/10 at 17:17:54", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "18", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp0"))
		insDoc(aux2, gLnk("R", "19", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp1"))
		insDoc(aux2, gLnk("R", "20", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp2"))
		insDoc(aux2, gLnk("R", "22", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp3"))
		insDoc(aux2, gLnk("R", "23", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp4"))
		insDoc(aux2, gLnk("R", "24", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp5"))
		insDoc(aux2, gLnk("R", "25", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp6"))
		insDoc(aux2, gLnk("R", "27", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_samp7"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "Aclarkii_Leyte_noAPR_Cf29,APY_65_2010-01-12.htm#12_01_10at17_17_54_comp_sum_numAll"))
