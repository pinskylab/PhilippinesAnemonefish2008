
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_2009-03-09 CebuLeyt noACH_C1.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 11/03/09 at 10:31:24", "Aclarkii_2009-03-09%20CebuLeyt%20noACH_C1.htm#11_03_09at10_31_24"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_2009-03-09%20CebuLeyt%20noACH_C1.htm#11_03_09at10_31_24_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_2009-03-09%20CebuLeyt%20noACH_C1.htm#11_03_09at10_31_24_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_2009-03-09%20CebuLeyt%20noACH_C1.htm#11_03_09at10_31_24_pairw_diff"))
