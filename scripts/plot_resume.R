library("argparse")
parser <- ArgumentParser(description='Plot resume log')
parser$add_argument('path', type="character",help="working directory path")

args <- parser$parse_args()

logfiles = list.files(path=args$path, pattern="*resume.log")
names = gsub(".resume.log", "", logfiles)

myList <- lapply(
    logfiles,
    function(x) {
    	print(x)
        temp = read.table(paste(args$path,x, sep="/"))
        temp$V3
        
    }
)
 
mytable = do.call(data.frame, myList)
names(mytable) = names

subtable = mytable[c(1,3) ,]

png("all.resume.png", width=1024, height=800)
barplot(as.matrix(subtable) / colSums(subtable))
dev.off()