distill <- function(file, gray = F){
  if(!gray){
    system(paste0("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.4 -sColorConversionStrategy=/CMYK -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=false -dAutoRotatePages=/None -sOutputFile=./temp.pdf ",file))
  }else{
    system(paste0("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.4 -sColorConversionStrategy=/Gray -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=false -dAutoRotatePages=/None -sOutputFile=./temp.pdf ",file))
  }
  system(paste0("rm ",file))
  system(paste0("mv ./temp.pdf ",file))
  if(!gray){
    system(paste0("convert -density 600 ",file," ",gsub(".pdf",".png",file)))
  }else{
    system(paste0("convert -colorspace gray -density 600 ",file," ",gsub(".pdf",".png",file)))
  }

}