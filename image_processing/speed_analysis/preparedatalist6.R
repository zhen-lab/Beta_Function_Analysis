# prep. .txt files named as genotype_number
# eg. wt_01.txt, wt_02.txt, ... e5_08.txt

#roi name for analysis
targetroiname<-"RIM"

#rois for anteroposterioraxis
#anterior
roiname1<-"nr"
#postrior
roiname2<-"RIM"

#colorrane for heatmap etc
#colorrange<-c(3,8)
#colorrange<-c(2,9)
#GCaMP3
colorrange<-c(0,0.5)

#range for trace figure
#tracerange<-c(3,10)
#GCaMP3
tracerange<-c(0,0.5)
#horizontalline<-5
horizontalline<-0.15
# 16x obj um/pix with  4x binning
scale<-0.39375*4
# 63x 0.1 um/pix? 4xbin
#scale<-0.1*4
# 40x 4xbin
#scale<-0.1575*4

#TRUE or FALSE
rotationapaxis<-FALSE
#rotationapaxis<-TRUE
 
#divideopposit <-FALSE
divideopposit<-TRUE



################################################################################## 
extractdata<- function(x, ...) {
	#prep sublist strings
	sublist<-NULL
	if(length(c(...))!=0)
	{
		sublist<-paste("$",paste(..., sep="$"),sep="")
	}
	#below line is old slow code
	#returnobjct<-NULL
	#this is faster
	returnobjct<-as.list(rep(NA,length(x)))
	for(i in 1:length(x))
	{
		#append is slow
		#returnobjct <-append(returnobjct,list(eval(parse(text=paste("x[[i]]",sublist,sep="")))))

		returnobjct[i] <-list(eval(parse(text=paste("x[[i]]",sublist,sep=""))))
		#x[[i]]$manualdata$crudegenotype
	}
	returnobjct
	#sublist
}


# This function take a vector and min max value of range.
# then return a vector consist with color code correspond with value.\
# out of ranged values are black.
getLUTvec <- function(x, min, max)
{
	#prepare color code
	startpoint<-2/3
	#red
	endpoint<-0
	hvec<-seq(startpoint, endpoint,len=255)
	temp<-1:255
	#function of (temp^3)/(255^3)/3 is decided so that looks nice for	 both normal and green-red blind. it will make darker col as near max
	colorcode<-hsv(h= hvec,s= 1,v= 1-(temp^3)/(255^3)/3)
	#
	colorrange<-c(min,max)
	rasio_range<-round((x-colorrange[1])/(colorrange[2]-	colorrange[1])*255)
	#black. out of range
	ratio_col<-rep("#000000",length(x))
	#white NA
	ratio_col[na.action(na.omit(x))]<-"#FFFFFF"
	for(i in 1:length(rasio_range))
	{
		if(!is.na(rasio_range[i]))
		{
			if(rasio_range[i]>0 && rasio_range[i]<256)
			{
				ratio_col[i]<-colorcode[rasio_range[i]]
			}
		}
	}
	return(ratio_col)
}


#################################
readDVTdata <-function(samplename, divideopposit=FALSE,...){	samplename<<-samplename		# to magage DVTracer10, read first 2lines and see if the column length is equal or not 	headers<-readLines(paste(samplename,".txt",sep=""),n=2)	roinamesstuff<-0	# if equal, it is older than v10	if(length(strsplit(headers[1],",")[[1]])==length(strsplit(headers[2],",")[[1]]))	{		dataframe<<-read.table(paste(samplename,".txt",sep=""),header=TRUE,sep=",")		header <-colnames(dataframe)		#the 1st column is time point		roinamesstuff<-1	}	# if not it is ver 10 or above	else	{		roiswidthheight<<-headers[1]		dataframe<<-read.table(paste(samplename,".txt",sep=""),header=TRUE,sep=",",skip=1)		header <-colnames(dataframe)		#time, x,y,zpos are exist before the rois		roinamesstuff<-4	}	roinames<-""	for(i in (1+roinamesstuff):(length(header)-1))	{		tempname<-strsplit(header[i],"_")[[1]][2]		frag=0		for(j in 1:length(roinames))		{			if(roinames[j]==tempname)	 		{	 			frag =1	 		}		}		if(frag==0)		{			if(tempname!="bg")			{				roinames<-cbind(roinames, tempname)			}		}	}	roinames<-roinames[2:length(roinames)]#eliminate the 1st factor, ""	roinames<<-roinames# left one is grobal	dm<<-matrix(0,nrow=nrow(dataframe),ncol=length(roinames))	for(i in 1:length(roinames))	{		dm[,i]<<-eval(parse(text=paste("(dataframe$l_", roinames[i],"_mean-dataframe$l_bg_mean)/(dataframe$r_", roinames[i],"_mean-dataframe$r_bg_mean)",sep="")))	}		if(divideopposit)	{		dm<<-1/dm;		}		colnames(dm)<-roinames	par(mar = c(2, 2, 0, 0), ps=8)	matplot(filter(dm,rep(1,5)/5),xaxs="i",ylab="",type="l",lty=1,...)	par(usr = c(0, 1, 0, 1),ps=12)	legend(0.8,0.9, roinames,bty='n',text.col=(1:length(roinames)))	legend(0.02,0.1, samplename,bty='n')	#dev.copy(png, file=paste(samplename,"ratiograph.png", sep=""),width=500, height=100)	#dev.off()}

############################# with each day data
# store data produced by DVT into list, and process it later

##########################################################
#get data files name
samplenameslist<-system("ls | grep '.txt'",intern = TRUE)
#eliminate memo.txt
for(i in 1:length(samplenameslist))
{
	if(samplenameslist[i]=="memo.txt")
	{
		samplenameslist[i]<-list(NULL)
	}
}
samplenamesvec<-unlist(samplenameslist)
samplenamesvec

##########################################################
#put each data into list
# This part can be used for common application
data_list<-NULL
for(i in 1:length(samplenamesvec))
{	
	samplename<-strsplit(samplenamesvec[i],".txt")[[1]]
	readDVTdata(samplename, divideopposit)
	#dev.copy(png, file=paste(samplename,"boxcarmat.png", sep=""),width=500, height=150)
	#dev.off()
	rawdatalist<-list(crudegenotype=strsplit(samplename,"_")[[1]][1],samplenumber=strsplit(samplename,"_")[[1]][2],DVTdata=dataframe,roinames=roinames,rawratiomatrix=dm, roiswidthheight = roiswidthheight)
	data_list<-append(data_list, list(list(rawdatalist=rawdatalist)))
}
save(data_list, file=paste(format(Sys.time(), "%Y%m%d_%H_%M"),"data_list.Rdata",sep=""))

##########################################################

#make a ratio matrix about the targetroi.

data_list2<-data_list
#targetroiname<-"AVA"

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list2))
{
	correctedmatrix<-extractdata(data_list2, "rawdatalist","rawratiomatrix")[[i]]
	nrowofeach<-cbind(nrowofeach ,nrow(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)

which(nrowofeach == maxrow)

aROIs<-matrix(NA, ncol=length(data_list2),nrow= maxrow)
correctedmatrix<-NULL
#for(i in 1:length(lapply(extractdata(data_list, "rawdatalist","roinames"),"==","AVE")))
for(i in 1:length(data_list2))
{

	correctedmatrix<-extractdata(data_list2, "rawdatalist","rawratiomatrix")[[i]]
	logicalvector<-unlist(lapply(extractdata(data_list2, "rawdatalist", "roinames"), "==", targetroiname)[[i]])
	#logicalvector<-unlist(lapply(extractdata(data_list2, "rawdatalist", "roinames"), "==", "AVAAVE")[[i]])
	#aROIs <-cbind(aROIs,correctedmatrix[,c(logicalvector)])
	aROIs[1:nrow(correctedmatrix),i] <- correctedmatrix[,c(logicalvector)]
}

##########################################################
#boxplot

cludegenotypevec<-unlist(extractdata(data_list2, "rawdatalist","crudegenotype"))
cludegenotypevec
typeslength<-nlevels(factor(cludegenotypevec))
typeslevels<-levels(factor(cludegenotypevec))
maxcasenum<-max(summary(factor(cludegenotypevec)))
meandataframe<-NULL
for(i in 1: typeslength)
{
	tempmeanvec<-apply(aROIs[,unlist(extractdata(data_list2, "rawdatalist","crudegenotype"))==typeslevels[i]],2,mean,na.rm=TRUE)
	meandataframe <-rbind(meandataframe, data.frame(mean= tempmeanvec, genotype=rep(typeslevels[i],length(tempmeanvec))))
		
}

par(bty="l")
plot(meandataframe$genotype, meandataframe$mean)
dev.copy(png, file=paste("meanboxplot.png", sep=""),width=400, height=480)
dev.off()


#n.s.
#t.test(meandataframe[meandataframe$genotype=="hpis190",]$mean, meandataframe[meandataframe$genotype=="hp428;hpls190",]$mean)

#n.s.
#t.test(meandataframe[meandataframe$genotype=="hpis190",]$mean, meandataframe[meandataframe$genotype=="nca(lf);hpls190",]$mean)


##########################################################
#heat map
#AVA D3
#colorrange<-c(3,8)

startpoint<-2/3
#red
endpoint<-0
hvec<-seq(startpoint, endpoint,len=255)
temp<-1:255
#function of (temp^3)/(255^3)/3 is decided so that looks nice for both normal and green-red blind
colorcode<-hsv(h= hvec,s= 1,v= 1-(temp^3)/(255^3)/3)

quartz()
#par(mfcol=c(typeslength,1),family="serif",mar=c(2,1,0,0))
par(mfcol=c(typeslength,1),mar=c(2,1,0,0))
#loop for genotype
for(i in 1: typeslength)
{
	tempbooleanvec<- unlist(extractdata(data_list2, "rawdatalist", "crudegenotype")) == typeslevels[i]
	emptymatrix<-matrix(NA, nrow=nrow(aROIs),ncol= maxcasenum)
	tempimagematrix<-filter(aROIs[, tempbooleanvec],rep(1,10)/10)
	emptymatrix[,1:ncol(tempimagematrix)]<-tempimagematrix
	imagematrix<-emptymatrix
	#par(mar = c(3, 2, 0, 2), ps=24)
	par(mar = c(1, 0, 0, 0), ps=24)
	image(0:nrow(imagematrix), 1:ncol(imagematrix), imagematrix, 	zlim= colorrange,col= colorcode, xlab="",ylab="", axes=FALSE)
}

dev.copy(png, file=paste("heatmap.png", sep=""),width=470, height=300)
dev.off()


#color bar
quartz()
colormatrix<-matrix((colorrange[1]*100):(colorrange[2]*100),nrow=1)
ticklabels <-(colorrange[2]-colorrange[1])*axTicks(4)+colorrange[1]
par(mar = c(2, 0, 2, 2))
image(colormatrix, zlim= colorrange*100,col= colorcode, xaxt="n", xlab="", ylab="")
axis(4, ticklabels, at= axTicks(4))
dev.copy(png, file=paste("colorcodebar.png", sep=""),width=40, height=400)
dev.off()


##########################################################


# show trace 
#tracerange<-c(3,8)
#horizontalline<-5
quartz()
par(mfcol=c(maxcasenum, typeslength),mar=c(2,1,0,0))
#loop for genotype
for(i in 1: typeslength)
{
	tempbooleanvec<- unlist(extractdata(data_list2, "rawdatalist", "crudegenotype")) == typeslevels[i]
	tempratiomatrix<-filter(aROIs[, tempbooleanvec],rep(1,10)/10)
	
	for(j in 1:maxcasenum)
	{
		par(mar = c(0, 1.5, 0, 0), ps=8,bty="l", mgp=c(1,0.5,0))
		if(j<=ncol(tempratiomatrix))
		{
			# what look like if caldium insensitive one has baseline around 6
			#plot(tempratiomatrix[,j]*1.5,type="l",col=1, ylim=c(5,11),xlab="",ylab="", axes=FALSE)
			plot(tempratiomatrix[,j],type="l",col=1, ylim=c(tracerange[1],tracerange[2]),xlab="",ylab="", axes=FALSE)
		axis(2)
		abline(h= horizontalline)
		}
		# beyond the case until max case
		else
		{
			#plot nothing
			plot(0,type="n",axes=FALSE)
		}
		par(usr = c(0, 1, 0, 1),ps=12)
		legend(0.0,0.55, j,bty='n')
	}
	par(usr = c(0, 1, 0, 1),ps=12)
	legend(0.1,0.55, typeslevels[i],bty='n')
	
}
dev.copy(png, file=paste("rawtrace.png", sep=""),width=800, height=1000)
dev.off()









##########################################################
# tracking data
#
#load("20100818data_list.Rdata")
#data_list2<-data_list

outputfigure<-TRUE
#roiname1<-"AVA"
#roiname2<-"c"

#return fperiod and rperiod
frperioddetection<-function(arg, threathold)
{
	#the argument is velocitywfr. raw or smoothen.
	#NOT have done. and threathold. lower than this treat as not moving
	
	velocitywfr <- arg	
	#fperiod plot
	#which return index. forward period
	# this filter(1,10) seems cause backperiod shift to left.
	# because backmotion is faster than forward, so at the bit eary poit velocity became 0 if running averaged velocity was used.
	#findex<-which(filter(velocitywfr,rep(1,10)/10)>0)
	
	#test firter(3)
	#findex<-which(filter(velocitywfr,rep(1,3)/3)>0)
	findex<-which(velocitywfr>0)
	
	#if diff(findex)>1, it mean change of dirrectin or end of time period
	endoff<-c(findex[1:(length(findex)-1)][diff(findex)>1], findex[length(findex)])
	startoff<-c(rev(findex)[1:(length(findex)-1)][diff(rev(findex))<(-1)], findex[1])
	startoff<-rev(startoff)
	#end processing
	if(endoff[1]< startoff[1])
	{
		startoff<-c(1,startoff)
	}
	if(endoff[length(endoff)]<startoff[length(startoff)])
	{
		endoff <-c(endoff,length(velocitywfr))
	}
	fperiod<-cbind(startoff, endoff)
	
	
	
	#rindex<-which(filter(velocitywfr,rep(1,3)/3)<0)
	rindex<-which(velocitywfr<0)
	
	#if diff(findex)>1, it mean change of dirrectin or end of time period
	endofr<-c(rindex[1:(length(rindex)-1)][diff(rindex)>1], rindex[length(rindex)])
	startofr<-c(rev(rindex)[1:(length(rindex)-1)][diff(rev(rindex))<(-1)], rindex[1])
	startofr<-rev(startofr)
	#end processing
	if(endofr[1]< startofr[1])
	{
		startofr<-c(1,startofr)
	}
	if(endofr[length(endofr)]<startofr[length(startofr)])
	{
		endofr <-c(endofr,length(velocitywfr))
	}
	rperiod<-cbind(startofr, endofr)


	

	list(fperiod, rperiod)	
}




#return absposx1, absposy1, absposx2, absposy2, rotangledir, velocitywfr
processvel<-function(testnum)
{
	workingDVTdata<-data_list2[[testnum]]$rawdatalist$DVTdata
	samplename <- paste(data_list2[[testnum]]$rawdatalist	$crudegenotype,data_list2[[testnum]]$rawdatalist$samplenumber,sep="_")
	dm<-data_list2[[testnum]]$rawdatalist$rawratiomatrix
	roinames<-data_list2[[testnum]]$rawdatalist$roinames
	roiswidthheight<- data_list2[[testnum]]$rawdatalist$roiswidthheight
	targetroiindex<-which(roinames==targetroiname)
	framepersec <- 1000/workingDVTdata$time[2]
	
	stagedata<-data.frame(x= workingDVTdata$xpos,y= workingDVTdata$ypos)
	stagedata$x<- -stagedata$x-min(-stagedata$x)
	stagedata$y<- -stagedata$y-min(-stagedata$y)
	# 63x 0.1 um/pix? 4xbin
	#scale<-0.1*4
	# 10x obj, 0.63?
	#scale<-0.63*4
	# 16x obj
	#scale<-0.39375*4
	# measuring showed actualy its about 0.47.? check if it change
	#scale<-0.47*4
	#seems not better than 0.39....
	# 40x obj
	#scale<-0.1575*4
	splitroiswidthheight<-unlist(strsplit(roiswidthheight,","))
	#110404 need to substitute some sign for Lin
	needtoescape<-grep("/",splitroiswidthheight)
	if(length(needtoescape)>0)
	{
		for(i in 1:length(needtoescape))
		{
			splitroiswidthheight[needtoescape[i]]<-gsub("/",".",splitroiswidthheight[needtoescape[i]])
		}
	}

	
	width1<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname1)+1])
	height1<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname1)+2])
	width2<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname2)+1])
	height2<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname2)+2])
	#cat("splitroiswidthheight",splitroiswidthheight,"\n")
	
	
	roi1xinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname1,"_x+ width1/2-168/2",sep="")))
	roi1yinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname1,"_y+ height1/2-256/2",sep="")))
	roi2xinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname2,"_x+ width2/2-168/2",sep="")))
	roi2yinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname2,"_y+ height2/2-256/2",sep="")))

	#use runmed to reduce wobbling 	
	# unit is micro meter	
	absposx1<-stagedata$x+ runmed(roi1xinfield,3)*scale
	absposy1<-stagedata$y+ runmed(roi1yinfield,3)*scale
	absposx2<-stagedata$x+ runmed(roi2xinfield,3)*scale
	absposy2<-stagedata$y+ runmed(roi2yinfield,3)*scale
	
	
	#if outputfigure==TRUE
	if(outputfigure)
	{
		ratiocolvec<-getLUTvec(dm[, targetroiindex], colorrange[1], colorrange[2])
		quartz()
		par(mar = c(2, 3, 0, 0), ps=8)
		#plot(absposx1, absposy1,type="o",pch=16,cex=0.5,col= 		ratiocolvec,xlim=c(0,500),ylim=c(500,0))
		
		
		plot(absposx1, absposy1,type="o",pch=16,cex=0.5,col= ratiocolvec,ylim=c(range(absposy1)[2],range(absposy1)[1]))
		points(absposx1[1], absposy1[1],pch="*",type="p",cex=3,col=1)
		points(absposx1[length(absposx1)], absposy1[length(absposx1)],pch="X",cex=1,type="p",col=1)

		dev.copy(png, file=paste(samplename,"posratio.png", sep=""),width=500, height=500)
		dev.off()
		dev.off()
	}

	
	
	
	
	
	
	
	
	
	#direcion and forward backward movement
	direction<-cbind(absposx1-absposx2,absposy1-absposy2)
	#plot(direction,xlim=c(-20,20),ylim=c(20,-20))
	if(outputfigure)
	{
	
		plot(direction, ylim=c(range(direction)[2], range(direction)[1]))
		dev.off()
	}
	angledir<-atan2(direction[,2], direction[,1])
	
	rotangledir<-angledir
	
	
	#rotation part
	if(rotationapaxis)
	{
		getkeyinput<-function(){
			temp<-readline(paste(samplename,"\n","Rotate roi2-roi1 vector? (input as clockwise degree)\n",sep=" "))
			while(is.na(as.numeric(temp)))
			{	
				temp <-readline(paste(samplename,"\n","Rotate roi2-roi1 vector? (input as clockwise degree)\n",sep=" "))
			}	
			temp
		}
	
		keyinput<-NA
		keyinput<-getkeyinput()
	
		dgreetheta<-as.numeric(keyinput)
		if(dgreetheta!=0)
		{
			theta<- dgreetheta/180*pi
			rotdir<-cbind(direction[,1]*cos(theta)-direction[,2]*sin(theta), direction[,1]*sin(theta)+direction[,2]*cos(theta))
			#plot(rotdir, ylim=c(range(rotdir)[2], range(rotdir)[1]))
	
			rotangledir <-atan2(rotdir[,2], rotdir[,1])
		}
	}#rotation part end
	
	#deal the apaxix failure
	# DIFFICALUT if bend tightly, worm change direction quickly.
	#but may help to check quality of apaxis detection
	#temp<-data_list2[[1]]$trackingdata$rotangledir
	#plot(temp,type="l")
	#diffofrad<-temp[2:length(rotangledir)]-temp[1:(length(rotangledir)-1)]
	diffofrad<-diff(rotangledir)
	#diffofrad[abs(diffofrad)>pi]<-abs(diffofrad[abs(diffofrad)>pi])-pi*2
	diffofradatan2<-atan2(sin(diffofrad),cos(diffofrad))
	#cat(length(diffofradatan2))
	#plot(diffofradatan2,type="l")
	#which(abs(diffofradatan2)>pi/2)
	if(outputfigure)
	{
		quartz()
		par(mfcol=c(2,1),mar = c(2, 3, 0, 0), mgp =c(2,1,0),ps=18)
		plot(rotangledir,type="o", ylab ="ap axis angle (radian)",ylim=c(-pi, pi))
		plot(diffofradatan2,type="o",ylab ="delta ap axis (radian)",ylim=c(-pi, pi))
		abline(h= pi/2,col=8)
		abline(h= -pi/2,col=8)
		dev.copy(png, file=paste(samplename,"apaxis.png", sep=""),width=1000, height=500)
		dev.off()		
		dev.off()		
	}	
	
	
	
	movedir<-atan2(diff(absposy1), diff(absposx1))
	# in case of large ref roi, it better?
	#movedir<-atan2(diff(filter(absposy2,rep(1,10)/10)), diff(filter(absposx2,rep(1,10)/10)))
	# is this way correct? +120 and -80 = +200 ok..
	frcos<-cos(movedir-rotangledir[2:length(rotangledir)])
	frfrag<-sign(frcos)
	#um/frame x10=um/sec
	velocitywfr<-frfrag*sqrt(diff(absposy1)^2+diff(absposx1)^2)
	#store non filtered data	
	#velocitywfr<-frfrag*sqrt(diff(filter(absposx2,rep(1,3)/3))^2+ 	diff(filter(absposy2,rep(1,3)/3))^2)
	
	
	
	
	
	
	
	
	
	if(outputfigure)
	{
		quartz()
		par(mar = c(2, 3, 0, 0), ps=8)
		#plot(filter(velocitywfr,rep(1,10)/10),type="l",xaxs="i")
		plot(filter(velocitywfr,rep(1,3)/3),type="l",xaxs="i",ylim=c(-10,10))
		#plot(filter(velocitywfr,rep(1,60)/60),type="l",xaxs="i")
		abline(h=0)
		dev.off()
		#dev.off()
	}
	
	if(outputfigure)
	{
	
		quartz()
		par(mar = c(3, 2, 0.5, 0), mgp=c(2,1,0), ps=8)
		plot(density(velocitywfr* framepersec,na.rm=TRUE,bw=0.2* framepersec), bty="l",main="density velocity",xlab="velocity (um/sec)")
		abline(v=0)
		#par(usr = c(0, 1, 0, 1),ps=26)
		#legend(0,1,mean(velocitywfr, na.rm=TRUE),bty="n")
		#roughly parcent 
		#forward
		legend("bottomright",paste(round(length(velocitywfr[velocitywfr>0])/(length(velocitywfr)),digit=3)),bty="n",cex=3)
		#stop 
		length(velocitywfr[velocitywfr==0])/(length(velocitywfr))
		#back 
		legend("bottomleft",paste(round(length(velocitywfr[velocitywfr<0])/(length(velocitywfr)),digit=3)),bty="n",cex=3)
		legend("topleft",samplename,bty="n")
		dev.copy(png, file=paste(samplename,"velositydensity.png", sep=""),width=400, height=400)
		dev.off()
		dev.off()
	}
	
	
	
	
	
	#get fperiod, rperiod
	#110531 fixed. this is not using NA at fisrt
	#templist <-frperioddetection(velocitywfr, 0)
		
	templist <-frperioddetection(c(NA,velocitywfr), 0)
	
	if(outputfigure)
	{
		quartz()
		par(mar = c(2, 3, 0, 0), ps=8,bty="l")
	
		plot(filter(velocitywfr,rep(1,10)/10),type="l",xaxs="i",ylim=c(-7,7))
		abline(h=0)
		#lines(filter(dm[,1],rep(1,10)/10)-5,col=2)
		#to precise correration, first frame is eliminated
	
		#lines((filter(dm[2:length(dm[, targetroiindex]),1],rep(1,10)/10)-horizontalline)/(tracerange[2]-tracerange[1])*10,col=2)
		lines((filter(dm[, targetroiindex],rep(1,10)/10)-horizontalline)/(tracerange[2]-tracerange[1])*10,col=2)
		for(i in 1:nrow(templist[[1]]))
		{
			lines(c(templist[[1]][i,]),c(rep(5,2)),col=2)
		}
	
	
		#rperiod plot
		#which return index. forward period
		#rindex<-which(filter(velocitywfr,rep(1,10)/10)<0)

		#plot(filter(velocitywfr,rep(1,10)/10),type="l",xaxs="i")
		#abline(h=0)
		#lines(filter(dm[,1],rep(1,10)/10)-5,col=2)
		for(i in 1:nrow(templist[[2]]))
		{
			lines(c(templist[[2]][i,]),c(rep(-5,2)),col=2)
		}
		legend("topleft", samplename,bty="n")
		dev.copy(png, file=paste(samplename, targetroiname,"velratioperiods.png", sep=""),width=800, height=300)
		dev.off()
		dev.off()
	}	
	
	
	

	#return absposx1, absposy1, absposx2, absposy2, rotangledir, velocitywfr, fperiod, rperiod
	#110531 Here the velocity has NA at first frame,
	# but in processvel frperioddetection function did'nt used NA added vector. fixed at above
	list(trackingdata= data.frame(cbind(absposx1, absposy1, absposx2, absposy2, rotangledir, velocity=c(NA,velocitywfr))), forwardperiod=templist[[1]], reverseperiod=templist[[2]])

}











# do all sample

preparevelocityall<-function()
{
	for(j in 1:length(samplenamesvec))
	{
		cat(samplenamesvec[j], "\n")
		#get absposx1, absposy1, absposx2, absposy2, rotangledir, velocitywfr
		templist<-processvel(j)
		data_list2[[j]]<-append(data_list2[[j]], templist)
		
		
		# for velocity csv export 
		samplename <- paste(data_list2[[j]]$rawdatalist	$crudegenotype,data_list2[[j]]$rawdatalist$samplenumber,sep="_")
		dm<-data_list2[[j]]$rawdatalist$rawratiomatrix
		roinames<-data_list2[[j]]$rawdatalist$roinames
		targetroiindex<-which(roinames==targetroiname)
		write.csv(data.frame(velocity=templist$trackingdata$velocity, ratio= dm[,targetroiindex]), file=paste(samplename,"velratio.csv",sep=""),row.names = FALSE)
		
		
		maxrow<-max(nrow(templist$forwardperiod), nrow(templist$reverseperiod))
		timepointsdata<-matrix(NA, ncol=4, nrow= maxrow)
		timepointsdata[1:nrow(templist$forwardperiod),1:2]<-templist$forwardperiod
		timepointsdata[1:nrow(templist$reverseperiod),3:4]<-templist$reverseperiod
		colnames(timepointsdata)<-c("fstart","fend","rstart","rend")
		
		write.csv(data.frame(timepointsdata), file=paste(samplename,"periodstimepoints.csv",sep=""),row.names = FALSE)

	}
	
	save(data_list2, file=paste(targetroiname,"data_list2.Rdata",sep=""))

}#end preparevelocity<-function()

#run above function
preparevelocityall()