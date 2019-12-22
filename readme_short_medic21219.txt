


Columns:

medication: the input data from questionnaires
MatchtoBNF2: the approximate matches to BNF inputs
All_Class : all corresponding medication class 
1-10:individual medication class
NumberMatchtoBNF2: the number of approximate matches of individual medication



Details
the data received 

character strings stripped from  the “tail end” of the input medication are in appendix A. These characters were nuanced to participants and likely to impede matching as they were not likely to be along the actual medication in the standard words or spellings such as in BNF


BNF product and chemical names columns formed the “seed” or “dictionary” of standard medication list. For the two columns, separately, duplicates were removed while keeping all corresponding subparagraphs (whether duplicated or not). To improve level of matching stringency and to account for usage of salt compound of some of the medications (with two words e.g.Atenolol Chloride ), the product names columns were splitted by whitespace, characters after space were removed , again while keeping the subparagraph column. All three “clean” (non-duplicated product name with subpara, non-duplicated chemical name with subpara, and non-duplicated space-splitted product name with subpara) 2-column sets were merged and again duplicates removed(not in the subparagraph)

BNF subparagraphs were used as medication class in output

For approximate matches, multiple algorithms and multiple scoring parameters and weights were trialled before using the current one in the scripts with distance varied generally between 70-90% similarity (actual value depending on algorithm being trialled,weights for 4 major ways of difference,  and length of the two strings to match e.g 0.15 for jato-winkler and 1.5 for levenshtein)