
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (AFG5-8,11-14,16,1,17-21 APR_Cf29 090208.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 08/02/09 at 22:47:59", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_47_59"))
	insDoc(aux1, gLnk("R", "Settings", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_47_59_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_47_59_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_47_59_amova"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_47_59_pairw_diff"))
	aux1 = insFld(foldersTree, gFld("Run of 08/02/09 at 22:53:25", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_53_25"))
	insDoc(aux1, gLnk("R", "Settings", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_53_25_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_53_25_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_53_25_amova"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "AFG5-8,11-14,16,1,17-21%20APR_Cf29%20090208.htm#08_02_09at22_53_25_pairw_diff"))
