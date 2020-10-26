
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Documents%20and%20Settings/Administrator/My%20Documents/PopGen%20Programs/Arlequin311/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_Cebu_noAC1359_2010-01-12.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 12/01/10 at 17:17:34", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "7", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp0"))
		insDoc(aux2, gLnk("R", "8", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp1"))
		insDoc(aux2, gLnk("R", "9", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp2"))
		insDoc(aux2, gLnk("R", "10", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp3"))
		insDoc(aux2, gLnk("R", "11", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp4"))
		insDoc(aux2, gLnk("R", "13", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp5"))
		insDoc(aux2, gLnk("R", "14", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp6"))
		insDoc(aux2, gLnk("R", "15", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp7"))
		insDoc(aux2, gLnk("R", "16", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp8"))
		insDoc(aux2, gLnk("R", "17", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_samp9"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "Aclarkii_Cebu_noAC1359_2010-01-12.htm#12_01_10at17_17_34_comp_sum_numAll"))