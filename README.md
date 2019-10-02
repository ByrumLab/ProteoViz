# ProteoViz
Interactive tool for Phosphoproteomics <br>


Set of R scripts for the analysis and visualization of Phosphoproteomics data sets. <br>

Script 1: Creates metadata that links with the MaxQuant database search results <br>
Script 2: Normalizes and runs limma differential expression on the protein lysate samples <br>
Script 3: Normalizes phosphosites to the protein expression and runs limma differential expression <br>
Script 4: Creates the input data file needed to run PTMsig analysis <br>
Script 5: Runs PTMsig/ssGSEA <br>
Script 6: Runs EGSEA gene set enrichment analysis <br>
Script 7A-7C: Creates the interactive Shiny dashboard for visualization

MaxQuant txt files are saved in the txt folder <br>
Scripts are in the src folder <br>
Create a Sample_metadata.tsv and contrast_matrix.tsv file associated with your experiment to run the code

# References:

1.	M. E. Ritchie, B. Phipson, D. Wu, Y. Hu, C. W. Law, W. Shi and G. K. Smyth, Nucleic Acids Research, 2015, 43, e47-e47. https://academic.oup.com/nar/article/43/7/e47/2414268

2.	K. Krug, P. Mertins, B. Zhang, P. Hornbeck, R. Raju, R. Ahmad, M. Szucs, F. Mundt, D. Forestier, J. Jane-Valbuena, H. Keshishian, M. A. Gillette, P. Tamayo, J. P. Mesirov, J. D. Jaffe, S. A. Carr and D. R. Mani, Mol Cell Proteomics, 2019, 18, 576-593. https://www.mcponline.org/content/18/3/576

3.	M. Alhamdoosh, C. W. Law, L. Tian, J. M. Sheridan, M. Ng and M. E. Ritchie, F1000Res, 2017, 6, 2010. https://f1000research.com/articles/6-2010
