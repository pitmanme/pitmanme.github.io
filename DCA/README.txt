DCA Tutorial in Development, by M Pitman - Version 1.0
———————————————————————————


——————————————————————
1. Software Required
——————————————————————
1.1- Matlab w/ Bioinformatics toolbox

1.2- HMMER (downloadable version is in the folder sent or visit http://hmmer.janelia.org) This maps the specific protein sequence of interest to the Hidden Markov Model (HMM) for the family of proteins in Pfam.

1.3- Chimera (Installer attached) Used to visualize the Direct Couplings.

1.4- mfDCA scripts (attached or can be obtained from http://dca.ucsd.edu/DCA/DCA.html)


———————————————————————————————————
2. Files Required Before Starting
———————————————————————————————————
All files used and generated are attached as well as examples from Faruck Morcos et al. However, it might be better to generate unique files to check for errors. 


2.1- FASTA file of the individual protein sequence. Go to PDB -> find protein (ex. 1AAR) -> Sequence tab on top -> download from FASTA link in center of page
/*do not close the webpage, more info on FASTA format in Note 4.1 */


2.2- Still in PDB, go to the Annotations tab and scroll to “Protein Family Annotation” -> obtain Pfam Accession code (PF00240 used for Ubiquitin)->

Go to http://pfam.xfam.org and in “JUMP TO” enter Accession code -> Go to Alignments tab on left -> Under Format an alignment select (Full or other of choice, FASTA format, Tree, Inserts lower case, Gaps as “.” or “-“ (mixed))

This will provide you with the Multiple Sequence Alignment (MSA) of the protein’s family. 
/*do not close webpage*/


2.3- In Pfam, go to the Curation & model tab on the left and then download the raw HMM for the family. For ubiquitin, the file I used is ubiquitin.hmm.txt


2.4- PDB of your structure


——————————————
3. PROCEDURE
——————————————
3.1- In Matlab set the path for  dca.m to perform mean-field DCA. Run dca and call according to the main function. In my example: 
	
	dca(‘MSA_ubiquiin.fasta’, ‘ubiquitin.DI’)

3.2- Sort the output file by the last column with a sequence distance greater than 4 (this may need to be changed)

	awk ‘$2-$1>4’ ubiquitin.DI|sort -g -k 4 -r > Ubiquitin_ranked.DI

3.3- Create mapping of specific protein onto family. On the PDB page where the Accession code for Pfam is found, a picture of this mapping is found. First create files accessed by HMMER
	
	hmmpress ubiquitin.hmm.txt

Then scan to create mapping with the FASTA file of the single protein

	hmmscan -o output_prefix —notextw HMMfile protein_FASTA_file
	hmmscan -o 1AAR_A_scan -notextw ubiquitin.hmm.txt 1AAR_A.fasta.txt

/* More information on the output file in Note 4.2 */

3.4- Combine output file with ranked DI pairs (see pg 64 of Faruck Morcos et al). This may be another function of hmmscan (unsure about this step)


3.5- Visualization. Create a 80x2 matrix variable for matlab. Mine is called ubiq_1AAR_top80.  This uses the top 80 contacts generated which I believe come from the output of Procedure step 2. Each row should be AA DI contact pairs. Also, it is stated the contacts should be mapped to the protein. To do this, I shifted the residue numbers by +5 based on the output of 3.3. 


3.6- Use plotDCAmap.m script with new variable. You must provide the path for heatmaps in the script (files attached). Call and Run with 

	[Mat] = plotDCAmap(ubiq_1AAR_top80,[], [1 76] , 1, 0)

I didn’t compare with a native structure yet. Just wanted to make sure it works. It generates a file attached but there are some errors flagged in matlab. 


3.7- Generate predicted contacts for Chimera. The dca_package provides GeneratePseudobonds.bash. The input of the top 80 is a txt file similar to the variable matrix created in Procedure step 5. There should be one space between each residue pair per row. 

	
	bash GeneratePseudobonds.bash offset Ubiquitin_contacts_top80.txt 1AAR_pairs 1AAR
	bash GeneratePseudobonds.bash offset top_contacts_matrix.txt filename_prefix chain_ID

It is probably better to use only the top30 by creating another 30x2 txt file. When I compared the example files provided, I think this is the right file to use for input in this step. 
This generates two pseudo bond files one *.dat and the other *.cmd. 


3.8- Visualize contacts in Chimera. In chimera fetch PDBID. Load Pseudobond*.dat file (or Tools -> Depiction -> Pseudobond reader -> select file). Run command pseudobond*.cmd (or File -> Open -> select file).

Edit attributes/add labels as desired.


——————————
4. NOTES
——————————

4.1- 

Example of FASTA format:

>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE

/* There should be a heading but some variation is acceptable (/ vs |) */


4.2-
Interpreting hmmscan outputs 
More info on:

	http://manpages.ubuntu.com/manpages/saucy/man1/hmmscan.1.html     
	http://hmmer.janelia.org/help/result#seqmatch

E values are significant at approx. < 0.001, but ideally should be << 0.001

“Conditional E-value - This is the E-value that the inclusion and reporting significant thresholds that are measured against (if defined as E-values). The conditional E-value is an attempt to measure the statistical significance of each domain, given that it has already been decided that the target sequence is a true homolog. It is the expected number of additional domains or hits that would be found with a domain/hit score this big in the set of sequences reported in the top hits list, if those sequences consisted only of random nonhomologous sequence outside the region that sufficed to define them as homologs.

Independent E-value - This is the significance of the sequence in the whole database search, if this were the only domain/hit that had been identified. If this E-value is not good, but the full sequence E-value is good, this is a potential red flag. Weak hits, none of which are good enough on their own, are summing up to lift the sequence up to a high score.”  (http://hmmer.janelia.org/help/result#seqmatch)

Bit score: The website does not specify the cut off for significance so I will search for that but bbased on examples ~30 (or greater) should be fine.

+ means that the residues are similar AA inserts

Example (expand to see format):

Alignments for each domain:
	  == domain 1  score: 122.2 bits;  conditional E-value: 2.8e-40
                      	                EETTSEEEEEEEETTSBHHHHHHHHHHHHTSTGGGEEEEETTEEETTTSBCGHHTTSTTEEEEEEESSS CS
 Query line,family ubiquitin-> 	      1 ktlsgktitlevepsdtvselKekieekegipvdqqrLiysGkvLeDdrtlsdyniqdgstlhlvlrlr 69
                                        ktl+gktitlevepsdt++++K+ki++kegip+dqqrLi++Gk+LeD+rtlsdyniq++stlhlvlrlr
Target-> 1AAR:A|PDBID|CHAIN|SEQUENCE  6 KTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLR 74
                                      8******************************************************************86 PP

4.3- For more information refer to “Direct Coupling Analysis for Protein Contact Prediction” Ch. 5, Faruck Morcos et al. (2014)


Okay, thanks for reading! ┏((＝￣(ｴ)￣=))┛




