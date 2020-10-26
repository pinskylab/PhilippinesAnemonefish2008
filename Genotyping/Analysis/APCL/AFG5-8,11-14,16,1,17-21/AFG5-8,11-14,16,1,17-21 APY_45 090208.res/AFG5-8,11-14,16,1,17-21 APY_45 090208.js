
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (AFG5-8,11-14,16,1,17-21 APY_45 090208.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 08/02/09 at 23:03:47", "AFG5-8,11-14,16,1,17-21%20APY_45%20090208.htm#08_02_09at23_03_47"))
	insDoc(aux1, gLnk("R", "Settings", "AFG5-8,11-14,16,1,17-21%20APY_45%20090208.htm#08_02_09at23_03_47_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "AFG5-8,11-14,16,1,17-21%20APY_45%20090208.htm#08_02_09at23_03_47_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "AFG5-8,11-14,16,1,17-21%20APY_45%20090208.htm#08_02_09at23_03_47_amova"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "AFG5-8,11-14,16,1,17-21%20APY_45%20090208.htm#08_02_09at23_03_47_pairw_diff"))
