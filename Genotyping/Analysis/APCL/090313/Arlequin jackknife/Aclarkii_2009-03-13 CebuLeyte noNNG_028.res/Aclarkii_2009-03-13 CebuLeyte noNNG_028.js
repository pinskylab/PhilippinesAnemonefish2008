
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_2009-03-13 CebuLeyte noNNG_028.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 15/03/09 at 17:15:22", "Aclarkii_2009-03-13%20CebuLeyte%20noNNG_028.htm#15_03_09at17_15_22"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_2009-03-13%20CebuLeyte%20noNNG_028.htm#15_03_09at17_15_22_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_2009-03-13%20CebuLeyte%20noNNG_028.htm#15_03_09at17_15_22_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_2009-03-13%20CebuLeyte%20noNNG_028.htm#15_03_09at17_15_22_pairw_diff"))
