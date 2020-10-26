setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/MattGribble")

afd2 = read.csv("../Extractions/AFD2 090622.csv", header=TRUE, row.names=1, na.string=c("Contam", "Null"))
colnames(afd2) = 1:12
afd2$row = row.names(afd2)
afd6 = read.csv("../Extractions/AFD6 090622.csv", header=TRUE, row.names=1, na.string="Null")
colnames(afd6) = 1:12
afd6$row = row.names(afd6)

todo = read.csv("MGsamplesToDo 090703.csv")

todo$num = gsub("PRBI", "", todo$Sample)

# turn afd2 and afd6 into long lists
cols = c("01", "02", "03","04","05","06","07","08","09","10","11","12")
a2 = data.frame(plate = "AFD2", col= "01", row=afd2$row, num=afd2[1:8,1])
for(i in 2:12){
	a2 = rbind(a2, data.frame(plate = "AFD2", col=cols[i], row=afd2$row, num=afd2[1:8,i]))
}
a6 = data.frame(plate="AFD6", col= "01", row=afd6$row, num=afd6[1:8,1])
for(i in 2:12){
	a6 = rbind(a6, data.frame(plate="AFD6", col=cols[i], row=afd6$row, num=afd6[1:8,i]))
}

# Match needed samples to location
todo = merge(todo, a2, all.x=T, by="num")
	todo$row = as.character(todo$row)
	todo$col = as.character(todo$col)
	todo$plate = as.character(todo$plate)
	todo$row[is.na(todo$row)] = ""
	todo$col[is.na(todo$col)] = ""
	todo$plate[is.na(todo$plate)] = ""
todo = merge(todo, a6, all.x=T, by="num")
	todo$row.y = as.character(todo$row.y)
	todo$col.y = as.character(todo$col.y)
	todo$plate.y = as.character(todo$plate.y)
	todo$row.y[is.na(todo$row.y)] = ""
	todo$col.y[is.na(todo$col.y)] = ""
	todo$plate.y[is.na(todo$plate.y)] = ""

todo$row = paste(todo$row.x, todo$row.y, sep="")
todo$col = paste(todo$col.x, todo$col.y, sep="")
todo$plate = paste(todo$plate.x, todo$plate.y, sep="")
dim(todo)

# Remove duplicates
todo = todo[!duplicated(paste(todo$Marker, todo$num)),]
dim(todo)
dim(todo)/96

i = order(todo$Marker, todo$col, todo$row)
todo = todo[i,c("Marker", "plate", "col", "row", "num", "Reason")]

table(todo$Marker)


# Setup the PCRs
pcrAC1578_1 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
cols = colnames(pcrAC1578_1)
rows = rownames(pcrAC1578_1)
i = which(todo$Marker=="AC1578" & todo$plate=="AFD2")
todo[i,]
for(j in i){
	pcrAC1578_1[rows==todo$row[j], cols==todo$col[j]] = paste(todo$num[j], "\n", todo$plate[j], " ", todo$col[j], todo$row[j], sep="")
}
write.csv(pcrAC1578_1, "AF16.M.csv")


pcrAC1578_2 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
cols = colnames(pcrAC1578_2)
rows = rownames(pcrAC1578_2)
i = which(todo$Marker=="AC1578" & todo$plate=="AFD6")
todo[i,]
for(j in i){
	pcrAC1578_2[rows==todo$row[j], cols==todo$col[j]] = paste(todo$num[j], "\n", todo$plate[j], " ", todo$col[j], todo$row[j], sep="")
}
write.csv(pcrAC1578_2, "AF17.M.csv")

pcrNNG_012_1 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
cols = colnames(pcrNNG_012_1)
rows = rownames(pcrNNG_012_1)
i = which(todo$Marker=="NNG_012" & todo$plate=="AFD2")
todo[i,]
for(j in i){
	pcrNNG_012_1[rows==todo$row[j], cols==todo$col[j]] = paste(todo$num[j], "\n", todo$plate[j], " ", todo$col[j], todo$row[j], sep="")
}
write.csv(pcrNNG_012_1, "AF18.M.csv")

pcrNNG_012_2 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
cols = colnames(pcrNNG_012_2)
rows = rownames(pcrNNG_012_2)
i = which(todo$Marker=="NNG_012" & todo$plate=="AFD6")
todo[i,]
for(j in i){
	pcrNNG_012_2[rows==todo$row[j], cols==todo$col[j]] = paste(todo$num[j], "\n", todo$plate[j], " ", todo$col[j], todo$row[j], sep="")
}
write.csv(pcrNNG_012_2, "AF19.M.csv")

pcr5 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="AC1359")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col=1
row=1
for(i in 1:split){
	pcr5[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}


i = which(todo$Marker=="APR_Cf29")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+4
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr5[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
write.csv(pcr5, "AF20.M.csv")
