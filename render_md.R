library(rmarkdown)
render("ENDOMETRIAL/endometrial.Rmd", md_document(variant = "gfm"),output_file = "README.md",output_dir = "ENDOMETRIAL")
render("BIRTHWEIGHT/birthweight.Rmd", md_document(variant = "gfm"), output_file = "README.md", output_dir = "BIRTHWEIGHT")
render("GRAPHS/graphs.Rmd", md_document(variant = "gfm"),output_file = "README.md",output_dir = "GRAPHS")
