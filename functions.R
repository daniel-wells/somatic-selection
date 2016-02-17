# If the file exists already, archive it by moving to archive folder with timestampd name so new file can be created in it's place
archive.file <- function(filename){
	if (file.exists(filename)){
		print(paste("Archiving old file:",filename))
	if (file.exists("archive")==FALSE){dir.create("archive", showWarnings = TRUE)}
	new.filename <- paste("archive/",format(file.mtime(filename),"%Y-%m-%d.%H-%M-%S."),filename,sep="")	
	file.rename(filename, new.filename)
	}else{
		print(paste("No file to archive:",filename))
	}
}