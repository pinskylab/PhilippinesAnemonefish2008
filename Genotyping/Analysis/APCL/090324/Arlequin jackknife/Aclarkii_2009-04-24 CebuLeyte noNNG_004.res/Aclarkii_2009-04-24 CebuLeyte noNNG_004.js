
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin%203.11/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Aclarkii_2009-04-24 CebuLeyte noNNG_004.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 24/04/09 at 10:13:50", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50"))
	insDoc(aux1, gLnk("R", "Settings", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_pairw_diff"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "7", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp0"))
		insDoc(aux2, gLnk("R", "8", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp1"))
		insDoc(aux2, gLnk("R", "9", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp2"))
		insDoc(aux2, gLnk("R", "10", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp3"))
		insDoc(aux2, gLnk("R", "11", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp4"))
		insDoc(aux2, gLnk("R", "13", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp5"))
		insDoc(aux2, gLnk("R", "14", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp6"))
		insDoc(aux2, gLnk("R", "15", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp7"))
		insDoc(aux2, gLnk("R", "16", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp8"))
		insDoc(aux2, gLnk("R", "17", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp9"))
		insDoc(aux2, gLnk("R", "18", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp10"))
		insDoc(aux2, gLnk("R", "19", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp11"))
		insDoc(aux2, gLnk("R", "20", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp12"))
		insDoc(aux2, gLnk("R", "22", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp13"))
		insDoc(aux2, gLnk("R", "23", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp14"))
		insDoc(aux2, gLnk("R", "24", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp15"))
		insDoc(aux2, gLnk("R", "25", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp16"))
		insDoc(aux2, gLnk("R", "27", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_samp17"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "Aclarkii_2009-04-24%20CebuLeyte%20noNNG_004.htm#24_04_09at10_13_50_comp_sum_numAll"))
