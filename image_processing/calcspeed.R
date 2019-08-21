# calculate a ROI speed. 
#20160815 Taizo

# 10x obj um/pix with  4x binning
#scale<-0.63*4
# 16x obj um/pix with  4x binning
# scale<-0.39375*4
# 63x 0.1 um/pix? 4xbin
scale<-0.1*4
# 40x 4xbin
# scale<-0.1575*4
# 20X 4bin
# scale<-0.315*4

outputfigure<-TRUE

#divideopposit <-FALSE
divideopposit<-TRUE


#################################
readDVTdata <-function(samplename, divideopposit=FALSE,...)
{
	samplename<<-samplename
	
	# to magage DVTracer10, read first 2lines and see if the column length is equal or not 
	headers<-readLines(paste(samplename,".txt",sep=""),n=2)
	roinamesstuff<-0
	# if equal, it is older than v10
	if(length(strsplit(headers[1],",")[[1]])==length(strsplit(headers[2],",")[[1]]))
	{
		dataframe<<-read.table(paste(samplename,".txt",sep=""),header=TRUE,sep=",")
		header <-colnames(dataframe)
		#the 1st column is time point
		roinamesstuff<-1
	}
	# if not it is ver 10 or above
	else
	{
		roiswidthheight<<-headers[1]
		dataframe<<-read.table(paste(samplename,".txt",sep=""),header=TRUE,sep=",",skip=1)
		header <-colnames(dataframe)
		#time, x,y,zpos are exist before the rois
		roinamesstuff<-4
	}
	roinames<-""
	for(i in (1+roinamesstuff):(length(header)-1))
	{
		tempname<-strsplit(header[i],"_")[[1]][2]
		frag=0
		for(j in 1:length(roinames))
		{
			if(roinames[j]==tempname)
	 		{
	 			frag =1
	 		}
		}
		if(frag==0)
		{
			if(tempname!="bg")
			{
				roinames<-cbind(roinames, tempname)
			}
		}
	}
	roinames<-roinames[2:length(roinames)]#eliminate the 1st factor, ""
	roinames<<-roinames# left one is grobal
	dm<<-matrix(0,nrow=nrow(dataframe),ncol=length(roinames))
	for(i in 1:length(roinames))
	{
		dm[,i]<<-eval(parse(text=paste("(dataframe$l_", roinames[i],"_mean-dataframe$l_bg_mean)/(dataframe$r_", roinames[i],"_mean-dataframe$r_bg_mean)",sep="")))
	}
	
	if(divideopposit)
	{
		dm<<-1/dm;	
	}
	
	colnames(dm)<-roinames

	par(mar = c(2, 2, 0, 0), ps=8)
	matplot(filter(dm,rep(1,5)/5),xaxs="i",ylab="",type="l",lty=1,...)
	par(usr = c(0, 1, 0, 1),ps=12)
	legend(0.8,0.9, roinames,bty='n',text.col=(1:length(roinames)))
	legend(0.02,0.1, samplename,bty='n')
	#dev.copy(png, file=paste(samplename,"ratiograph.png", sep=""),width=500, height=100)
	#dev.off()

}

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
#save(data_list, file=paste(format(Sys.time(), "%Y%m%d_%H_%M"),"data_list.Rdata",sep=""))


##########################################################
#return absposx1, absposy1, absposx2, absposy2, rotangledir, velocitywfr
processspeed<-function(testnum)
{
	workingDVTdata<-data_list[[testnum]]$rawdatalist$DVTdata
	samplename <- paste(data_list[[testnum]]$rawdatalist	$crudegenotype,data_list[[testnum]]$rawdatalist$samplenumber,sep="_")
	dm<-data_list[[testnum]]$rawdatalist$rawratiomatrix
	roinames<-data_list[[testnum]]$rawdatalist$roinames
	roiswidthheight<- data_list[[testnum]]$rawdatalist$roiswidthheight
	#targetroiindex<-which(roinames==targetroiname)
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

	
	#width1<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname1)+1])
	#height1<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname1)+2])
	#width2<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname2)+1])
	#height2<-as.numeric(splitroiswidthheight[which(splitroiswidthheight== roiname2)+2])
	width<-as.numeric(splitroiswidthheight[2])
	height<-as.numeric(splitroiswidthheight[3])
	#cat("splitroiswidthheight",splitroiswidthheight,"\n")
	
	
	#roi1xinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname1,"_x+ width1/2-168/2",sep="")))
	#roi1yinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname1,"_y+ height1/2-256/2",sep="")))
	#roi2xinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname2,"_x+ width2/2-168/2",sep="")))
	#roi2yinfield<-eval(parse(text=paste("workingDVTdata$l_", roiname2,"_y+ height2/2-256/2",sep="")))
	roixinfield<-eval(parse(text=paste("workingDVTdata$l_", roinames[1],"_x+ width/2-168/2",sep="")))
	roiyinfield<-eval(parse(text=paste("workingDVTdata$l_", roinames[1],"_y+ height/2-256/2",sep="")))

	#use runmed to reduce wobbling 	
	# unit is micro meter	
	#absposx1<-stagedata$x+ runmed(roi1xinfield,3)*scale
	#absposy1<-stagedata$y+ runmed(roi1yinfield,3)*scale
	#absposx2<-stagedata$x+ runmed(roi2xinfield,3)*scale
	#absposy2<-stagedata$y+ runmed(roi2yinfield,3)*scale
	absposx<-stagedata$x+ runmed(roixinfield,3)*scale
	absposy<-stagedata$y+ runmed(roiyinfield,3)*scale
	
	
	
	
	
	#um/frame x10=um/sec
	roispeed<-sqrt(diff(absposy)^2+diff(absposx)^2)
}
	

# do all sample

preparespeedall<-function()
{
	for(j in 1:length(samplenamesvec))
	{
		cat(samplenamesvec[j], "\n")
		#get absposx1, absposy1, absposx2, absposy2, rotangledir, velocitywfr
		tempspeed<-c(NA,processspeed(j))
		# for speed csv export 
		samplename <- paste(data_list[[j]]$rawdatalist	$crudegenotype,data_list[[j]]$rawdatalist$samplenumber,sep="_")
		dm<-data_list[[j]]$rawdatalist$rawratiomatrix
		roinames<-data_list[[j]]$rawdatalist$roinames
		#targetroiindex<-which(roinames==targetroiname)
		targetroiindex<-1
		if(outputfigure)
		{
			quartz()
			par(mar = c(0, 2, 0.5, 0.5), mgp=c(1,0.5,0), ps=8,bty="n",mfrow=c(2,1))
			ratiovec<-filter(dm[, targetroiindex],rep(1,3)/3)
			#ratiorange<-range(ratiovec,na.rm=T)
			plot(filter(dm[, targetroiindex],rep(1,3)/3),type="l",xaxs="i", col=2, axes=F,ylab="ratio")
			axis(2)
			speedvec<-filter(tempspeed,rep(1,3)/3)
			speedrange<-range(speedvec,na.rm=T)
			par(mar = c(2, 2, 0, 0.5),bty="l")
			plot(speedvec,type="l",xaxs="i",ylim=c(0, speedrange[2]),xlab="frame",ylab="speed (um/frame)")
			#abline(h=0)
			dev.copy(png, file=paste(samplename,"ratiospeed.png", sep=""),width=500, height=500)
			dev.off()
			dev.off()
		}
		write.csv(data.frame(speed= tempspeed, ratio= dm[,targetroiindex]), file=paste(samplename,"speedratio.csv",sep=""),row.names = FALSE)
	}
}

#run above function
preparespeedall()

