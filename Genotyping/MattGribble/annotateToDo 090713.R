setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/MattGribble")

afd2 = read.csv("../Extractions/AFD2 090622.csv", header=TRUE, row.names=1, na.string=c("Contam", "Null"))
colnames(afd2) = 1:12
afd2$row = row.names(afd2)
afd6 = read.csv("../Extractions/AFD6 090622.csv", header=TRUE, row.names=1, na.string="Null")
colnames(afd6) = 1:12
afd6$row = row.names(afd6)

todo = read.csv("MGsamplesToDo 090714.csv")

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
length(unique(todo$Marker))

i = order(todo$Marker, todo$col, todo$row)
todo = todo[i,c("Marker", "plate", "col", "row", "num", "Reason")]

table(todo$Marker)


# Setup the PCRs

pcr1 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="ACH_A11")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col=1
row=1
for(i in 1:split){
	pcr1[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
i = which(todo$Marker=="ACH_A3")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+12
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr1[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
write.csv(pcr1, "AF22.M.csv")

pcr2 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="ACH_A4")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col=1
row=1
for(i in 1:split){
	pcr2[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
i = which(todo$Marker=="ACH_A7")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+5
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr2[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
write.csv(pcr2, "AF24.M.csv")

pcr3 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="ACH_B9")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col=1
row=1
for(i in 1:split){
	pcr3[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
i = which(todo$Marker=="NNG_028")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+5
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr3[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
write.csv(pcr3, "AF25.M.csv")

pcr4 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="ACH_A8")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col=1
row=1
for(i in 1:split){
	pcr4[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
i = which(todo$Marker=="NNG_004")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+11
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr4[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
i = which(todo$Marker=="NNG_007")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+5
if(col >12){
	col=1
	row=row+1
}
for(i in 1:split){
	pcr4[row,col] = paste(now$num[i], "\n", now$plate[i], " ", now$row[i], now$col[i], sep="")
	col = col+1
	if(col >12){
		col=1
		row=row+1
	}
}
write.csv(pcr4, "AF26.M.csv")


pcr5 = matrix(nrow=8, ncol=12, dimnames = list(c("a","b","c","d","e","f","g","h"), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")), byrow=TRUE)
i = which(todo$Marker=="ACH_C1")
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
i = which(todo$Marker=="ACH_D1")
now = todo[i,]
now = now[order(now$plate, now$row, now$col),]
split = length(now$Marker)
col = col+12
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
write.csv(pcr5, "AF27.M.csv")

