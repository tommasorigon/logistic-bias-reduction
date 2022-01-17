library(rmarkdown)
render("ENDOMETRIAL/endometrial.Rmd", md_document(variant = "gfm"))
system("cp ENDOMETRIAL/endometrial.md ENDOMETRIAL/README.md")

render("BIRTHWEIGHT/birthweight.Rmd", md_document(variant = "gfm"))
system("cp BIRTHWEIGHT/birthweight.Rmd BIRTHWEIGHT/README.md")
