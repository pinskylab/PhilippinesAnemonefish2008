
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_sizes_2009-03-16.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 16/03/09 at 13:22:41", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_amova"))
		insDoc(aux2, gLnk("R", "FIS per pop", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_amova_AMOVA_FIS"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_pairw_diff"))
		insDoc(aux2, gLnk("R", "Locus by locus AMOVA", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41Loc_by_Loc_AMOVA"))
		insDoc(aux2, gLnk("R", "F-stat bootstraps", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_comp_sum_bootstrap"))
		insDoc(aux2, gLnk("R", "FIS per pop per locus", "Aclarkii_sizes_2009-03-16.htm#16_03_09at13_22_41_comp_sum_LBL_AMOVA_FIS"))
