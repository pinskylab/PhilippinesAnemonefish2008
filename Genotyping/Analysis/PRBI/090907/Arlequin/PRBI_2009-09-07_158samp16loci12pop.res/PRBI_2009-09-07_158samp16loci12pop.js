
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Documents%20and%20Settings/Lab/Desktop/Genetic%20Analysis%20Programs%20-%20Used/Arlequin%203.1/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (PRBI_2009-09-07_158samp16loci12pop.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 18/09/09 at 08:37:00", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00"))
	insDoc(aux1, gLnk("R", "Settings", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_amova"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "1", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp0"))
		insDoc(aux2, gLnk("R", "2", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp1"))
		insDoc(aux2, gLnk("R", "7", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp2"))
		insDoc(aux2, gLnk("R", "8", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp3"))
		insDoc(aux2, gLnk("R", "9", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp4"))
		insDoc(aux2, gLnk("R", "10", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp5"))
		insDoc(aux2, gLnk("R", "11", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp6"))
		insDoc(aux2, gLnk("R", "13", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp7"))
		insDoc(aux2, gLnk("R", "14", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp8"))
		insDoc(aux2, gLnk("R", "15", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp9"))
		insDoc(aux2, gLnk("R", "19", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp10"))
		insDoc(aux2, gLnk("R", "22", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_samp11"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_37_00_comp_sum_thetaH"))
	aux1 = insFld(foldersTree, gFld("Run of 18/09/09 at 08:41:33", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_41_33"))
	insDoc(aux1, gLnk("R", "Settings", "PRBI_2009-09-07_158samp16loci12pop.htm#18_09_09at08_41_33_run_information"))
