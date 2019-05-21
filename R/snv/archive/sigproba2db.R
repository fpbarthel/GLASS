

con <- DBI::dbConnect(odbc::odbc(), "GLASSv2")

tmp <- read.delim('Downloads/signatures_probabilities.txt', as.is = TRUE)

df <- tmp %>% gather(-(1:3),key="signature",value="proba") %>% filter(grepl("^Signature",signature))

df2 <- df %>% transmute(signature = as.numeric(gsub("Signature.","",signature)),
                        ref_context = Trinucleotide,
                        alt = substring(Substitution.Type,3,3),
                        substitution_type = Substitution.Type,
                        proba)

dbWriteTable(con, Id(schema="ref",table="signature_proba"), df2, overwrite = TRUE)
