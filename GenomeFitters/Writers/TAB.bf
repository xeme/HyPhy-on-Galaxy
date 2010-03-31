/* 	Written by SL Kosakovsky Pond in June 2007.	spond@ucsd.eduTAB.bf generates a tab-separated text file, with oneline of output per gene; invalid genes will be marked withID, user error (provided by the reader module). *//* these variables will be used by all three functions */_outputFilePath = "";_returnHeaders	= {}; /* stores the associative array of valid result keys,					     provided by a call to returnResultHeaders *//*---------------------------------------------------------*/function _prepareFileOutput (_outPath){	_outputFilePath = _outPath;		/* retrieve the descriptions of all results	   from the analysis module; store them for later */	   	_returnHeaders 	= returnResultHeaders(0);		/* write out the headers */	fprintf (_outputFilePath, CLEAR_FILE, KEEP_OPEN, _returnHeaders[0]);	for (_biterator = 1; _biterator < Abs(_returnHeaders); _biterator = _biterator + 1)	{		fprintf (_outputFilePath,"\t",_returnHeaders[_biterator]);	}	fprintf (_outputFilePath,"\n");	return 0;}/*---------------------------------------------------------*/function _processAGene (valid, _geneID, _errorTag){	if (valid)	{		returnValue = runAGeneFit (_geneID);		fprintf (_outputFilePath, returnValue[_returnHeaders[0]]);		for (_biterator = 1; _biterator < Abs(_returnHeaders); _biterator = _biterator + 1)		{			fprintf (_outputFilePath,"\t",returnValue[_returnHeaders[_biterator]]);		}		fprintf (_outputFilePath, "\n");	}	else	{		fprintf (_outputFilePath, _geneID,"\t", _errorTag, "\n");	}	_currentState = 0;	return 0;}/*---------------------------------------------------------*/function _finishFileOutput (dummy){	fprintf (_outputFilePath,CLOSE_FILE);	return 0;}