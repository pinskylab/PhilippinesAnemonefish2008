
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_GenAlEx_2009-03-09 onepop.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 10/03/09 at 15:43:32", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "All", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_samp0"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_43_32_comp_sum_numAll"))
	aux1 = insFld(foldersTree, gFld("Run of 10/03/09 at 15:44:09", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "All", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_samp0"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "Aclarkii_GenAlEx_2009-03-09%20onepop.htm#10_03_09at15_44_09_comp_sum_numAll"))
