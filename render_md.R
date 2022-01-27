library(rmarkdown)
render("ENDOMETRIAL/endometrial.Rmd", md_document(variant = "gfm"),output_file = "README.md",output_dir = "ENDOMETRIAL")
render("BIRTHWEIGHT/birthweight.Rmd", md_document(variant = "gfm"), output_file = "README.md", output_dir = "BIRTHWEIGHT")
render("HIGH-DIMENSIONAL-SYNTHETIC/simulation-studies.Rmd", md_document(variant = "gfm"), output_file = "README.md", output_dir = "HIGH-DIMENSIONAL-SYNTHETIC/")
