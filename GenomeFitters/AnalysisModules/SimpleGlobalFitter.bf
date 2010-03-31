/* 
	Written by SL Kosakovsky Pond in June 2007.
	spond@ucsd.edu

SimpleGlobalFitter.bf takes a single coding aligment, fits a 
MG94x(used define nucleotide biases model) to the alignment, 
and returns the following quantities per gene:

	[NUMBER]"GENE" 				: the index of the gene
	[NUMBER]"BP"   				: length of the alignment (in bp)
	[NUMBER]"S_sites"			: (expected) # of synonymous sites in the alignment
	[NUMBER]"NS_sites"			: (expected) # of non-synonymous sites in the alignment
	[NUMBER]"Stop_codons"		: # of alignment columns which contain at least one 
								  stop codon
	[NUMBER]"LogL"				: Log likelihood score of the codon model
	[NUMBER]"omega"				: Global estimate of dN/dS for the gene
	[STRING]"omega_range"		: The 95% lower and upper confidence bounds 
								  on omega, estimated by profile likelihood
								  The format is lower_bound-upper_bound
								  (e.g. "0.8-1.1")
	
	[NUMBER]"AC"				: A<->C substitution rate (relative to the A<->G rate)
	[NUMBER]"AT"				: A<->T substitution rate (relative to the A<->G rate)
	[NUMBER]"CG"				: C<->G substitution rate (relative to the A<->G rate)
	[NUMBER]"CT"				: C<->T substitution rate (relative to the A<->G rate)
	[NUMBER]"GT"				: G<->T substitution rate (relative to the A<->G rate)
	
	[STRING]"Tree"				: The Newick string with the tree for this alignment
	
	
Inputs (when the first analysis is executed):
	
	1). [STRING] 	Genetic code: the name of the genetic code (valid choices are
					listed below; they are case-sensitive!)
					
			Universal
			Vertebrate mtDNA
			Yeast mtDNA
			Mold/Protozoan mtDNA
			Invertebrate mtDNA
			Ciliate Nuclear
			Echinoderm mtDNA
			Euplotid Nuclear
			Alt. Yeast Nuclear
			Ascidian mtDNA
			Flatworm mtDNA
			Blepharisma Nuclear
				
	2). [STRING] 	A six-character encoding of the nucleotide substitution model
	(e.g http://www.hyphy.org/cgi-bin/yabb/YaBB.pl?board=DM1;action=display;num=1106770520;start=0#0)
	
Dependancies (in Utility):
	
		chooseGeneticCode.def
		dSdNTreeTools.ibf
		CodonTools.bf
		GrabBag.bf
		MG94custom.mdl
		NJ.bf (if treeString is empty)
		
External variables dependenacies
		
		[STRING] 	treeString: stores the tree string to be used by the analysis
		[DATASET]	ds: the coding alignment to be processed

*/

/*---------------------------------------------------------*/

VERBOSITY_LEVEL			   = -1;
COUNT_GAPS_IN_FREQUENCIES  = 0;
_totalFilesProcessed 	   = 0;

/*---------------------------------------------------------*/

function returnResultHeaders (dummy)
{
	_analysisHeaders = {};
	_analysisHeaders[0]  = "GENE";
	_analysisHeaders[1]  = "BP";
	_analysisHeaders[2]  = "S_sites";
	_analysisHeaders[3]  = "NS_sites";
	_analysisHeaders[4]  = "Stop_codons";
	_analysisHeaders[5]  = "LogL";
	_analysisHeaders[6]  = "omega";
	_analysisHeaders[7]  = "omega_range";
	_analysisHeaders[8]  = "AC";
	_analysisHeaders[9]  = "AT";
	_analysisHeaders[10] = "CG";
	_analysisHeaders[11] = "CT";
	_analysisHeaders[12] = "GT";
	_analysisHeaders[13] = "Tree";

	return _analysisHeaders;
}

/*---------------------------------------------------------*/

function runAGeneFit (myID)
{
	fprintf (stdout, "[SimpleGlobalFitter.bf on GENE ", myID, "]\n");
	taxonNameMap = {};
	
	for (k=0; k<ds.species; k=k+1)
	{
		GetString 		(thisName, ds,k);
		shortName 		= (thisName^{{"\\..+",""}})&&1;
		taxonNameMap[shortName] = thisName;
		SetParameter (ds,k,shortName);
	}

	DataSetFilter filteredData = CreateFilter (ds,1);
	_nucSites	 			   = filteredData.sites;
	
	if (_totalFilesProcessed == 0)
	/* first call to this function; include the NJ tree module */
	{
		/* read in utility files; see documentation for what the files do
		   and the source for all available batch funtions */
		
			/* this file will prompt for the genetic code, reading the input
			   from the overloaded standard input in the configuration file 
			   It populates the 64 state vector mapping (_Genetic_Code) codon to an 
			   (internal) aminoacid code, and a string of comma separated stop codons 
			   (e.g. "TAA,TAG,TGA" for the universal code) stored in GeneticCodeExclusions
			   
			*/
		ExecuteAFile 			("../Utility/chooseGeneticCode.def");
			   
			/* permits the construction of dS and dN trees */
		ExecuteAFile 			("../Utility/dSdNTreeTools.ibf");
			/* defines code needed to compute synonymous and non-synonymous
			   sites per codon */
		ExecuteAFile 			("../Utility/CodonTools.bf");
			   
		ExecuteAFile 			("../Utility/GrabBag.bf");

	}
	if (Abs(treeString))
	{
		givenTreeString = treeString;
	}
	else
	{
		if (_totalFilesProcessed==0)
		{
			/* read in the 6-character model specification for 
			   nucleotide biases */
			/* modelSpecString is the name of the variable used by
			   model defintion files */

			ExecuteAFile 			("../Utility/NJ.bf");
		}
		givenTreeString = InferTreeTopology (0);
		treeString		= "";
	}

	DataSetFilter filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

	if (_totalFilesProcessed==0)
	/* first call to this function; read in model description and 
	   execute the file which defines the substitution model */
	{
		fscanf					(stdin, "String",  modelSpecString);
		
		/* provide the necessary options for the model definition via 
		   overloaded standard input */
		   
		_MG94stdinOverload = {};
		_MG94stdinOverload ["0"] = "Global";
		_MG94stdinOverload ["1"] = modelSpecString;
		
		/* Execute the model definition file */
		ExecuteAFile 			("../Utility/MG94custom.mdl", _MG94stdinOverload);
	}
	else
	{
		/* these commands are executed in MG94custom.mdl when it is first 
		   included, but for subsequent runs, they must be done explicitly */
		   
		HarvestFrequencies 			  (observedFreq,filteredData,3,1,1);
		MULTIPLY_BY_FREQS 		    = PopulateModelMatrix ("MG94custom", observedFreq);
		vectorOfFrequencies 	    = BuildCodonFrequencies (observedFreq);
		Model MG94customModel 		= (MG94custom,vectorOfFrequencies,0);

	}

	Tree	codonTree		    = givenTreeString;
	LikelihoodFunction lf     = (filteredData,codonTree);

	Optimize 					(res,lf);

		/* 
		  this function is defined in CodonTools.bf and computes an associative array of total, synonymous
		   and non-synonymous sites in the DataSetFilter 
		*/
		
	_snsAVL	   = 				_computeSNSSites ("filteredData", _Genetic_Code, vectorOfFrequencies, 0);
	
		/* this function is defined in dSdNTreeTools.bf and computes an associative array of "codonTree"
		   scaled on various quantities, such as the expected number of substitution/codon, dS, dN etc
		   The current implementation does not actually use the result...
		   
		*/
	_cL		   =  				ReturnVectorsOfCodonLengths (ComputeScalingStencils (0), "codonTree");


	_returnMe = {};
	_returnMe ["GENE"]  			= myID;
	_returnMe ["LogL"]  			= res[1][0];
	_returnMe ["BP"] 				= _snsAVL ["Sites"];
	_returnMe ["S_sites"] 			= _snsAVL ["SSites"];
	_returnMe ["NS_sites"] 			= _snsAVL ["NSSites"];
	_returnMe ["Stop_codons"] 		= (_nucSites-filteredData.sites*3)$3;
	_returnMe ["AC"] 				= AC;
	_returnMe ["AT"] 				= AT;
	_returnMe ["CG"] 				= CG;
	_returnMe ["CT"] 				= CT;
	_returnMe ["GT"] 				= GT;
	
	/* MG94custom.mdl with the 'Global setting' defines 'R' to be the global dN/dS ratio */
	_returnMe ["omega"] 			= R;
	COVARIANCE_PARAMETER 			= "R";
	COVARIANCE_PRECISION 			= 0.95;
	CovarianceMatrix 				(cmx,lf);
	_returnMe ["omega_range"] 		= ""+cmx[0]+"-"+cmx[2];
	_returnMe ["Tree"] 				= Format(codonTree,0,1);

	_totalFilesProcessed = _totalFilesProcessed + 1;
	
	return _returnMe;
}

