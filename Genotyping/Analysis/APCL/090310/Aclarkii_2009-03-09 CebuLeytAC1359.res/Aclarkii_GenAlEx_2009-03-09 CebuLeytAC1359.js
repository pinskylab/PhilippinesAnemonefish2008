
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_GenAlEx_2009-03-09 CebuLeytAC1359.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 10/03/09 at 16:23:22", "Aclarkii_GenAlEx_2009-03-09%20CebuLeytAC1359.htm#10_03_09at16_23_22"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_GenAlEx_2009-03-09%20CebuLeytAC1359.htm#10_03_09at16_23_22_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_GenAlEx_2009-03-09%20CebuLeytAC1359.htm#10_03_09at16_23_22_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_GenAlEx_2009-03-09%20CebuLeytAC1359.htm#10_03_09at16_23_22_pairw_diff"))
