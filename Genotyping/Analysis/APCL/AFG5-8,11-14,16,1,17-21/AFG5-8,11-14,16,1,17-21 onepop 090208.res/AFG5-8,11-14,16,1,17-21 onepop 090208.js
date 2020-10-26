
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (AFG5-8,11-14,16,1,17-21 onepop 090208.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 08/02/09 at 22:38:27", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27"))
	insDoc(aux1, gLnk("R", "Settings", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "All", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_samp0"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_numAll"))
		insDoc(aux2, gLnk("R", "Allelic range", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_allRange"))
		insDoc(aux2, gLnk("R", "Garza-Williamson index", "AFG5-8,11-14,16,1,17-21%20onepop%20090208.htm#08_02_09at22_38_27_comp_sum_GW"))
